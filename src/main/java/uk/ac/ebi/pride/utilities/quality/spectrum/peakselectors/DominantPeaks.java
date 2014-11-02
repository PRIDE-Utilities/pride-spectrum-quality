/*
 * Created on Feb 2, 2004
 */
package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;


import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.spectrum.common.SpectrumUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * PeakSelector that returns "dominant" peaks: every peak
 * that's above 20% of the average of the 5 largest peaks
 * is returned.
 *
 * @author ypriverol
 */
public class DominantPeaks implements PeakSelector{

    @Override
	public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException {

        Spectrum analyzeSpectrum = (Spectrum) spectrum.clone();

        analyzeSpectrum = SpectrumUtils.sortSpectrumByIntensity(analyzeSpectrum, false);

		int n = 0;
		double sum = 0;
		// sum up to 5 biggest ones
		for(int i = 0; i < analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray().length && i < 5;i++)		{
			sum += (analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray()[(spectrum.getIntensityBinaryDataArray().getDoubleArray().length-1 - i)]);
			n++;
		}

		if(n == 0){
			return null;
		}

		double criterium = 0.2 * sum/n;

		// apply cutoff
        List<Double> mzList  = new ArrayList<Double>();
        List<Double> intList = new ArrayList<Double>();

        for(int i = 0; i < analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray().length; i++){
            if(analyzeSpectrum.getIntensityBinaryDataArray().getDoubleArray()[i] >= criterium){
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
		return "Dominant peaks: bigger than 20% of the average of the 5 biggest peaks";
	}

    @Override
    public String getLabel() {
        return "Dom. peaks";
    }
}
