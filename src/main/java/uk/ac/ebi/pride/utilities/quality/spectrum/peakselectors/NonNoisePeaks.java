package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import java.text.NumberFormat;
import java.util.ArrayList;

/**
 * A PeakSelector that returns peaks that are probably not
 * noise. This implementation returns a peak when it's above
 * the mean plus a certain configurable multiple of the
 * standard deviation of all peaks.
 * SPECTRUM NOT PARTITIONED
 *
 * @author M. Vogelzang
 */
public class NonNoisePeaks implements PeakSelector{

	protected double zLimit;

	protected int deletionLimit;

	public NonNoisePeaks(double zLimit, int deletionLimit){
		this.zLimit = zLimit;
		this.deletionLimit = deletionLimit;
	}

	public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException {

        ArrayList numbers = new ArrayList();
		//Spectrum analyzedSpectrum = (Spectrum) spectrum.clone();
        Spectrum analyzedSpectrum = new Spectrum(spectrum);

		double[][] peakList = analyzedSpectrum.getMassIntensityMap();

		for (int i = 0; i < peakList.length; i++){
			numbers.add(new Float(Math.log(peakList[i][1] + 1E-8)));
		}

		RobustMath.NormalDistributionParameters ndp = RobustMath.robustEstimate((Float[])numbers.toArray(new Float[0]), deletionLimit, zLimit);

		numbers.clear();
		for (int i = 0;i <peakList.length;i++)
			if(Math.log(peakList[i][1] + 1E-8) > ndp.mean + zLimit*ndp.stddev)
			{
				numbers.add(new Integer(i));
			}

		double[][] newPeakList = new double[numbers.size()][];

		for(int i = 0;i<numbers.size();i++)
			newPeakList[i] = peakList[((Integer)numbers.get(i)).intValue()];
		    analyzedSpectrum.setMassIntensityMap(newPeakList);

		return analyzedSpectrum;
	}

    @Override
	public String getDescription(){
        NumberFormat nf = NumberFormat.getInstance();
		return "Peaks determined to be 'real' because their intensity is bigger than the mean + " + nf.format(zLimit) + " SD";
	}

    @Override
	public String getLabel(){
		return "Sig. peaks";
	}
}
