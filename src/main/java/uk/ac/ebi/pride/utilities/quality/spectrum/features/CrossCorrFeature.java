package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import java.util.HashMap;
import java.util.Map;

/**
 * CrossCorrFeature is the class that provides a method to calculate
 * an identification-independent crosscorrelation score. This
 * reflects how many pairs of peaks could be matched b- and y-
 * ions.
 *
 * originally developed by M. Vogelzang
 * @author ypriverol
 *
 */
public class CrossCorrFeature implements FeatureCalculator{

    private Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

	protected final static int TOLERANCE = 2;

	protected static int getSingleHits(double target, boolean[] peakThere, int start){
		int hits = 0;
		for(int i = 0;i<peakThere.length;i++)
			if(peakThere[i]){
			boolean matched = false;
			int mz = i+start;
			int j = (int)(target - TOLERANCE - mz - start + 0.5);
			if(j < 0) j =0;
			int endj = (int)(target + TOLERANCE - mz - start + 0.5);
			if(endj >= peakThere.length) endj = peakThere.length-1;
			for(;j<=endj;j++)
				if(peakThere[j])
					hits++;
		}
		return hits;
	}

	protected static int getDoubleHits(double target, boolean[] peakThere, int start){

        int hits = 0;
		for(int i = 0;i<peakThere.length;i++)
			if(peakThere[i])
			{
				boolean matched = false;
				double mass = (i+start)-1;
				double secondMz = ((target * 2 - 2) - mass) + 1;
				int j = (int)(secondMz - TOLERANCE - start + 0.5); // what are the 0.5?
				if(j < 0) j =0;
				int endj = (int)(secondMz + TOLERANCE - start + 0.5);
				if(endj >= peakThere.length) endj = peakThere.length-1;
				for(;j<=endj;j++)
					if(peakThere[j])
						matched = true;
				if(matched)
				{
					hits ++;
					i++;
				}
			}
		return hits;
	}

	protected static int getTripleHits(double target, boolean[] peakThere, int start)
	{
		int hits = 0;
		for(int i = 0;i<peakThere.length;i++)
			if(peakThere[i])
			{
				boolean matched = false;
				int mz = i+start;
				double secondMass = (target * 3 - 3) - (mz-1);
				int j = (int)((secondMass+2)/2.0 - TOLERANCE - start + 0.5);
				if(j < 0) j =0;
				int endj = (int)((secondMass+2)/2.0 + TOLERANCE - start + 0.5);
				if(endj >= peakThere.length) endj = peakThere.length-1;
				for(;j<=endj;j++)
					if(peakThere[j])
						matched = true;
					if(matched){
						i++;
						hits++;
					}

			}
		return hits;
	}

    private static CrossCorrFeature instance = null;

    protected CrossCorrFeature(){
        features.put(SpectrumFeatureType.QUALSCORE_CROSSCORR_B_Y_IONS, 0);
    }

    public static CrossCorrFeature getInstance(){
        if(instance == null)
            instance = new CrossCorrFeature();
        return instance;
    }

    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

        double[][] peakList = spectrum.getMassIntensityMap();

        if (peakList.length == 0)
            return features;

        int minMz = (int) peakList[0][0];
        int maxMz = (int) (peakList[peakList.length - 1][0] + 1);

        boolean[] peakThere = new boolean[maxMz - minMz + 1];

        for (int i = 0; i < peakList.length; i++)
            peakThere[(int) (peakList[i][0] - minMz)] = true;

        int[] dev = new int[]{-10, -15, -20, -25, -30, 10, 15, 20, 25, 30};

        double precursorMZ = DataAccessUtilities.getPrecursorMz(spectrum);

        if (charge == 1) {

            int realHits = getSingleHits(precursorMZ, peakThere, minMz);

            Integer[] values = new Integer[dev.length];
            for (int i = 0; i < dev.length; i++)
                values[i] = new Integer(getSingleHits(precursorMZ + dev[i], peakThere, minMz));

            RobustMath.NormalDistributionParameters ndp = RobustMath.normalEstimate(values);

            double z;

            if (ndp.stddev < 0.001) {
                if (realHits < ndp.mean - 0.5)
                    z = -2;
                else if (realHits > ndp.mean + 0.5)
                    z = 2;
                else
                    z = 0;
            } else
                z = (realHits - ndp.mean) / ndp.stddev;

            features.put(SpectrumFeatureType.QUALSCORE_CROSSCORR_B_Y_IONS, z);

            return features;
        }

        // Check doubly / triply charged
        int realHits = getDoubleHits(precursorMZ, peakThere, minMz);

        Integer[] values = new Integer[dev.length];

        for (int i = 0; i < dev.length; i++)
            values[i] = new Integer(getDoubleHits(precursorMZ + dev[i], peakThere, minMz));

        RobustMath.NormalDistributionParameters ndp = RobustMath.normalEstimate(values);

        double z2;

        if (realHits == 0 && ndp.mean < 0.01)
            z2 = 0;
        else {
            if (ndp.stddev > 0.3)
                z2 = (realHits - ndp.mean) / ndp.stddev;
            else
                z2 = (realHits - ndp.mean) * 3;

            realHits = getTripleHits(precursorMZ, peakThere, minMz);
            for (int i = 0; i < dev.length; i++)
                values[i] = new Integer(getTripleHits(precursorMZ + dev[i], peakThere, minMz));

            ndp = RobustMath.normalEstimate(values);

            double z3;
            if (realHits == 0 && ndp.mean < 0.01)
                z3 = 0;
            else {
                if (ndp.stddev > 0.3)
                    z3 = (realHits - ndp.mean) / ndp.stddev;
                else
                    z3 = (realHits - ndp.mean) * 3;
            }

            if (z3 > z2) {
                features.put(SpectrumFeatureType.QUALSCORE_CROSSCORR_B_Y_IONS, z3);
                return features;
            } else {
                features.put(SpectrumFeatureType.QUALSCORE_CROSSCORR_B_Y_IONS, z3);
                return features;
            }
        }
        return features;
    }

}
