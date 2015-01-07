package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.ProcessingType;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import javax.swing.text.NumberFormatter;
import java.text.NumberFormat;
import java.util.ArrayList;

/**
 * Another implementation of PeakSelector that is meant to
 * return the peaks that are not noise. This implementation
 * does the same as NonNoisePeaks, except that it partitions
 * the spectrum into a number of m/z divisions of which
 * it separately takes peaks above a certain intensity percentile,
 * i.e. the peaks are sorted by intensity and e.g. the highest
 * 10% of all peaks are declared to be signal peaks.
 * @author M. Vogelzang
 */
public class NonNoisePeaks3 implements PeakSelector{

	protected double zLimit, percentileCutoff;

	protected int deletionLimit, maxPeaksPer1000Da;

	protected int divisions;

	public NonNoisePeaks3(int divisions, double zCutoff, int maxPeaks){
		this.divisions = divisions;  // the spectrum is divided into m/z-intervals since the intensity varies a lot between the center and the periphery of the spectry
		this.percentileCutoff = zCutoff;
		this.maxPeaksPer1000Da = maxPeaks;
	}

	public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException {

		ArrayList numbers = new ArrayList(), includeNumbers = new ArrayList();

        //Spectrum analyzeSpectrum = (Spectrum) spectrum.clone();
        Spectrum analyzeSpectrum = new Spectrum(spectrum);

		double[][] peakList = analyzeSpectrum.getMassIntensityMap();

		if (peakList.length <= divisions){
			return null;
		}

		int currentIndex = 0;

		double minMz = peakList[0][0], maxMz = peakList[peakList.length - 1][0];

		for (int k = 0; k < divisions; k++)
		{
			double currentMax = minMz + (maxMz - minMz) * (k + 1) / divisions;
			numbers.clear();
			int oldIndex = currentIndex;
			for (;
				currentIndex < peakList.length
					&& peakList[currentIndex][0] <= currentMax;
				currentIndex++)
			{
				numbers.add(new Float(Math.log(peakList[currentIndex][1] + 1E-8)));
			}
			RobustMath.NormalDistributionParameters ndp =
				RobustMath.estimatePercentiles(
					(Float[]) numbers.toArray(new Float[0]), percentileCutoff, maxPeaksPer1000Da * peakList.length / 1000);

            for (int i = oldIndex; i < currentIndex; i++)
				if (Math.log(peakList[i][1] + 1E-8)
					> ndp.mean)
				{
					includeNumbers.add(new Integer(i));
				}
		}

		double[][] newPeakList = new double[includeNumbers.size()][];

		for (int i = 0; i < includeNumbers.size(); i++)
			newPeakList[i] = peakList[((Integer) includeNumbers.get(i)).intValue()];
		     analyzeSpectrum.setMassIntensityMap(newPeakList);

		return analyzeSpectrum;
	}
    @Override
	public String getDescription(){
        return ProcessingType.NONNOSE3_PEAK.getTitle();
	}

    @Override
	public String getCode(){
        return ProcessingType.NONNOSE3_PEAK.getCode();
	}
}
