package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.ProcessingType;

/**
 * The simplest PeakSelector: all peaks that are input are outputted
 * too.
 * 
 * @author ypriverol
 */
public class AllPeaks implements PeakSelector{

    public Spectrum transform(Spectrum spectrum){
		return spectrum;
	}

    @Override
	public String getDescription(){
		return ProcessingType.ALL_PEAKS.getTitle();
	}

    @Override
    public String getCode() {
        return ProcessingType.ALL_PEAKS.getCode();
    }
}
