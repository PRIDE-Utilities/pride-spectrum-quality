package uk.ac.ebi.pride.utilities.quality.utils;

import uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors.*;

/**
 * All the processing methods we already implemented in the library.
 * @ypriverol
 */
public enum ProcessingType {

    ALL_PEAKS         ("ALLPEAK",        "Include all the peaks in the analysis", AllPeaks.class),
    DOM_PEAKS         ("DOMPEAK",        "Dominant peaks: bigger than 20% of the average of the 5 biggest peaks", DominantPeaks.class),
    DOM_DA_PEAKS      ("DOMDAPEAK",      "Dominant peaks: bigger than 20% of the average of the 5 biggest peaks", DominantPeaksPerDalton.class),
    ISO_PEAKS         ("ISOPEAK",        "Isotope signal, peaks which fit a mass-dependent isotope pattern predicted by a Poisson model", IsotopePeaks.class),
    NOISE_PEAKS       ("NOISEPEAK",      "Peaks determined to be noise because their intensity is smaller than the mean", NoisePeaks.class),
    NOISEINT_PEAK     ("NOISEINTPEAK",   "Noise peaks (Spectrum divided into divisions. Noise peaks are below intensity percentile", NoisePeaksByInterval.class),
    NONNOISE_PEAK     ("NONOISEPEAK",    "Peaks determined to be 'real' because their intensity is bigger than the mean", NonNoisePeaks.class),
    NONNOSE3_PEAK     ("NONOISE3PEAK",   "Signal peaks (Spectrum divided into divisions. Signal peaks are above intensity percentile: ", NonNoisePeaks3.class);

    private String code;
    private String title;
    private Class type;

    private ProcessingType(String code, String title, Class type) {

        this.code = code;
        this.title = title;
        this.type = type;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public Class getType() {
        return type;
    }

    public void setType(Class type) {
        this.type = type;
    }

    public String getCode() {
        return code;
    }

    public void setCode(String code) {
        this.code = code;
    }

}
