package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.spectrum.common.SpectrumUtils;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import java.util.ArrayList;
import java.util.List;

/**
 * A PeakSelector that intends to return only the peaks that are noise.
 * These are peaks that lower than the (robust) mean plus a certain
 * multiple of the standard deviation of all peaks.
 *
 * @author ypriverol
 */
public class NoisePeaks implements PeakSelector{

	protected double zLimit; // has a twofold function: for robust mean and as peak selection cutoff

	protected int deletionLimit;

	public NoisePeaks(double zLimit, int deletionLimit){
		this.zLimit = zLimit;
		this.deletionLimit = deletionLimit;
	}

	public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException {

        Spectrum analyzeSpectrum = new Spectrum(spectrum);

		analyzeSpectrum = SpectrumUtils.sortSpectrumByIntensity(analyzeSpectrum,false);
        List<Float> numbers = new ArrayList<Float>();

		for (int i = 0; i < analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray().length; i++){
			numbers.add(new Float(Math.log(analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray()[i] + 1E-8)));
		}

		RobustMath.NormalDistributionParameters ndp = RobustMath.robustEstimate(numbers.toArray(new Float[0]), deletionLimit, zLimit);

        // apply cutoff
        List<Double> mzList  = new ArrayList<Double>();
        List<Double> intList = new ArrayList<Double>();

        for (int i = 0;i < analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray().length;i++){
            if(Math.log(analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray()[1] + 1E-8) <= ndp.mean + zLimit*ndp.stddev)			{
                mzList.add(analyzeSpectrum.getMzBinaryDataArray().getDoubleArray()[i]);
                intList.add(analyzeSpectrum.getMzBinaryDataArray().getDoubleArray()[i]);
            }
        }

        double[] mzArr  = new double[mzList.size()];
        double[] intArr = new double[intList.size()];

        for(int i = 0; i < mzList.size();i++){
            mzArr[i]  = mzList.get(i);
            intArr[i] = intList.get(i);
        }

        analyzeSpectrum.getIntensityBinaryDataArray().setDoubleArray(intArr);
        analyzeSpectrum.getMzBinaryDataArray().setDoubleArray(mzArr);

		return analyzeSpectrum;
	}

    @Override
	public String getDescription(){
		return "Peaks determined to be noise because their intensity is smaller than the mean + " + zLimit;
	}

	public String getLabel(){
		return "Noise";
	}
}
