
package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import java.util.HashMap;
import java.util.Map;

/**
 * ParameterGen generates features, and is called by MasterFeatureGenerator
 * on different sets of peaks.
 * The features ParameterGen currently calculates are (please keep this documentation
 * in sync with the code at all times to avoid confusion):
 * <UL>
 * <LI>Number of peaks
 * <LI>Average intensity (log is taken to make the distribution more uniform)
 * <LI>Intensity standard deviation (log is taken here too)
 * <LI>"Normalized m/z range covered": the smallest range that includes 95% of all
 * intensities, divided by the precursorMz
 * <LI>"Normalized 50% m/z range covered": the same as above, but with 50% instead of
 * 95 %
 * <LI>Mean number of peaks / m/z: the total number of peaks divided by the 95%
 * m/z range covered. (log is taken)
 * <LI>Standard deviation of m/z gaps between consecutive peaks (log is taken)
 * <LI>Average number of peaks within distance &lt; 2 m/z units of another peak
 * </UL>
 *
 * originally developed by M. Vogelzang
 * @author ypriverol
 */
public class GeneralFeatures implements FeatureCalculator{

    private static GeneralFeatures instance;

    private Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    protected  GeneralFeatures(){
        features.put(SpectrumFeatureType.QUALSCORE_NUM_PEAKS,         0);
        features.put(SpectrumFeatureType.QUALSCORE_AVG_BY_INTENSITY,  0);
        features.put(SpectrumFeatureType.QUALSCORE_STD_INTENSITY,     0);
        features.put(SpectrumFeatureType.QUALSCORE_MZ_95_INTENSITY,   0);
        features.put(SpectrumFeatureType.QUALSCORE_MZ_50_INTENSITY,   0);
        features.put(SpectrumFeatureType.QUALSCORE_TIC_MZ,            0);
        features.put(SpectrumFeatureType.QUALSCORE_MASS_GAP,          0);
        features.put(SpectrumFeatureType.QUALSCORE_NEIGHGOR_2DA,      0);
    }

    public static GeneralFeatures getInstance(){
       if(instance == null)
           instance = new GeneralFeatures();
        return instance;
    }

    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

        if (spectrum == null || spectrum.getMassIntensityMap().length == 0)
            return features;

        double[][] peakList = spectrum.getMassIntensityMap();

        //Peak count
        int peakCount = peakList.length;

        Double[] peaks = new Double[peakList.length];

        // TIC total ion current
        double totalIntensity = 0;
        for (int i = 0; i < peakList.length; i++){
            peaks[i] = new Double(peakList[i][1]); // save intensities
            totalIntensity += peakList[i][1];
        }

        // normal distribution parameters
        RobustMath.NormalDistributionParameters ndp;
        ndp = RobustMath.normalEstimate(peaks);

        double avgIntensity      = ndp.mean;
        double sdIntensity       = ndp.stddev;
        double skewnessIntensity = ndp.skewness;

        // signal range
        double mzRange1 = getSmallestMassRangeContainingIntensity(peakList, 0.95 * totalIntensity);

        double mzRange2 = getSmallestMassRangeContainingIntensity(peakList, 0.50 * totalIntensity);

        // peak density
        double peakPerMz;
        double ticPerMz;

        if (mzRange1 <= 0.001)
            peakPerMz = 1;
        else
            peakPerMz = peakCount / mzRange1;

        // TIC per Da
        if (mzRange1 <= 0.001)
            ticPerMz = 1;
        else
            ticPerMz = totalIntensity / mzRange1;

        // analyse gap width between peaks
        peaks = new Double[peakList.length - 1];

        for (int i = 0; i < peakList.length - 1; i++)
            peaks[i] = new Double(peakList[i + 1][0] - peakList[i][0]);

        RobustMath.NormalDistributionParameters ndp2;
        ndp2 = RobustMath.normalEstimate(peaks);

        double sdMassGap;
        if (peakList.length == 1)
            sdMassGap = 0;
        else
            sdMassGap = ndp2.stddev;

        double outp_avgWithin2 = 0;

        /*
         *   average number of neighbours surrounding each peak within a mass range of +/- 2 Da
         *     illustration below
         *                   i   i  I     i  ii i  spectrum
         *                   i  Xi  I   X i  ii i  X: 2 Da limits  I: current peak  i: other peaks
         *                   1   b  e   4 5  67 8  peak numeration
         *
         *     e-b = 3-2 = 1  -> it means that the current peak I has one neighbour within +/- 2 Da
         */

        int begin = 0;
        int end = 0;
        for (int i = 0; i < peakList.length; i++){

            while (end < peakList.length - 1 && peakList[end + 1][0] - peakList[i][0] < 2) // if next peak is closer than 2 Da from current peak
                end++;
            while (begin < i && peakList[i][0] - peakList[begin][0] > 2) // if last peak is further than 2 Da from current peak
                begin++;

            outp_avgWithin2 += end - begin;  // [FR] I removed the plus 1, results in linear offset, no effect on LDA performance
        }

        outp_avgWithin2 /= peakList.length;

        double precursorMZ = DataAccessUtilities.getPrecursorMz(spectrum);

        features.put(SpectrumFeatureType.QUALSCORE_NUM_PEAKS, Math.sqrt(peakCount));
        features.put(SpectrumFeatureType.QUALSCORE_AVG_BY_INTENSITY,  Math.log(avgIntensity + 1));
        features.put(SpectrumFeatureType.QUALSCORE_STD_INTENSITY, Math.log(sdIntensity + 1));
        features.put(SpectrumFeatureType.QUALSCORE_MZ_95_INTENSITY, mzRange1 / precursorMZ);
        features.put(SpectrumFeatureType.QUALSCORE_MZ_50_INTENSITY, mzRange2 / precursorMZ);
        features.put(SpectrumFeatureType.QUALSCORE_TIC_MZ, Math.log(ticPerMz + 1E-4));
        features.put(SpectrumFeatureType.QUALSCORE_MASS_GAP, Math.log(sdMassGap + 1E-4));
        features.put(SpectrumFeatureType.QUALSCORE_NEIGHGOR_2DA, outp_avgWithin2);

        return features;


    }

	protected double[] getDefaults() {
		// defaults:returned when there is no data
		// should not be too much of an outlier
		return new double[] {0, 10, 10, 0, 0, -4, 4, 1 };
	}

	private double getSmallestMassRangeContainingIntensity(double[][] peakList, double targetIntensity) {

		double currentIntensitySum = peakList[0][1];
		int firstPeak = 0;
		int lastPeak = 0;

		while (currentIntensitySum < targetIntensity) {
                    currentIntensitySum += peakList[++lastPeak][1];
		}

		double smallestMassRange = peakList[lastPeak][0] - peakList[firstPeak][0];
                
		// smallest mass range in which 95 % of intensity is.
		while (lastPeak < peakList.length - 1) {
                    currentIntensitySum += peakList[++lastPeak][1];
                    while (firstPeak < peakList.length && (currentIntensitySum - peakList[firstPeak][1] >= targetIntensity)) {
                        currentIntensitySum -= peakList[firstPeak++][1];
                    }
                    if (peakList[lastPeak][0] - peakList[firstPeak][0] < smallestMassRange) {
                        smallestMassRange = peakList[lastPeak][0] - peakList[firstPeak][0];
                    }
		}
		return smallestMassRange;
	}

}
