package uk.ac.ebi.pride.utilities.quality.utils;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;

/**
 * @author ypriverol
 */
public enum SpectrumFeatureType {

    QUALSCORE_AA_MASS_DIFF                   ("QS-AAMAS-DIF",    "AA mass differences, number of", "Plain number of peak-distances that could be amino acids (normalized by dividing this number by the square of the total number of peaks)",Double.class),
    QUALSCORE_AA_MASS_ABUNDANCE_WEIGHTED     ("QS-AAMAS-AWT",    "AA mass differences, abundance-weighted",  "Feature 4 is roughly equivalent to feature 1 (number of AA distances) but it is an abundance-weighted version (FR)", Double.class),
    QUALSCORE_AA_MASS_TAG_AVG_LONGER         ("QS-AAMAS-TAGAVG", "Sequence tag, average longest",      "(Average) Based on chains of consecutive amino-acid-like distances: for every peak <I>i</I> this <sub>i</sub> will represent the maximum number of steps of amino-acid distance can be taken towards lower m/z values",Double.class),
    QUALSCORE_AA_MASS_TAG_LONGER             ("QS-AAMAS-TAGLNG", "Sequence tag, longest", "(Longest) Based on chains of consecutive amino-acid-like distances: for every peak <I>i</I> this <sub>i</sub> will represent the maximum number of steps of amino-acid distance can be taken towards lower m/z values",Double.class),

    QUALSCORE_Complement_B_Y_IONS_3vs3       ("QS-COMPL-3VS3",   "Parent mass tolerance 0.3 versus 3.0", "Parent mass tolerance 0.3 versus 3.0", Double.class),
    QUALSCORE_Complement_B_Y_IONS_1vs3       ("QS-COMPL-1VS3",   "Parent mass tolerance 1.0 versus 3.0", "Parent mass tolerance 1.0 versus 3.0", Double.class),
    QUALSCORE_Complement_B_Y_IONS_iso1       ("QS-COMPL-ISO1",   "Isotopes/background signal, parent_mass+1/parent_mass-10", "Isotopes/background signal, parent_mass+1/parent_mass-10", Double.class),
    QUALSCORE_Complement_B_Y_IONS_sgn        ("QS-COMPL-SIGN",   "Complements4 (main signal/background, PM+0/PM-10)", "Complements4 (main signal/background, PM+0/PM-10)", Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_SC ("QS-COMPL-3CHARGE-SC","Complements: simple count","Complements: simple count",Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_AB ("QS-COMPL-3CHARGE-AB","Complements: divide by averaged background","Complements: divide by averaged background",Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_SA ("QS-COMPL-3CHARGE-SA","Complements: subtract averaged background","Complements: subtract averaged background",Double.class),

    QUALSCORE_CROSSCORR_B_Y_IONS             ("QS-CROSS-CORR", "An indepent identification-independent crosscorrelation score", "An indepent identification-independent crosscorrelation score", Double.class),

    QUALSCORE_AMONIA_17_SC                   ("QS-AMONIA-17-SC","Neutral losses: ammonia (-17) simple count","Neutral losses: ammonia (-17) simple count",Double.class),
    QUALSCORE_AMONIA_17_AVGB                 ("QS-AMONIA-17-AVGB","Neutral losses: ammonia (-17) divided by averaged background","Neutral losses: ammonia (-17) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_17_AVGMB                ("QS-AMONIA-17-AVGMB","Neutral losses: ammonia (-17) count minus averaged background","Neutral losses: ammonia (-17) count minus averaged background",Double.class),
    QUALSCORE_AMONIA_18_SC                   ("QS-AMONIA-17-SC","Neutral losses: water (-18) simple count","Neutral losses: water (-18) simple count",Double.class),
    QUALSCORE_AMONIA_18_AVGB                 ("QS-AMONIA-18-AVGB","Neutral losses: ammonia (-18) divided by averaged background","Neutral losses: ammonia (-18) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_18_AVGMB                ("QS-AMONIA-18-AVGMB","Neutral losses: ammonia (-18) count minus averaged background","Neutral losses: ammonia (-18) count minus averaged background",Double.class),
    QUALSCORE_AMONIA_28_SC                   ("QS-AMONIA-17-SC","Neutral losses: water (-28) simple count","Neutral losses: water (-28) simple count",Double.class),
    QUALSCORE_AMONIA_28_AVGB                 ("QS-AMONIA-18-AVGB","Neutral losses: ammonia (-28) divided by averaged background","Neutral losses: ammonia (-28) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_28_AVGMB                ("QS-AMONIA-18-AVGMB","Neutral losses: ammonia (-28) count minus averaged background","Neutral losses: ammonia (-28) count minus averaged background",Double.class),

    QUALSCORE_ISOTOPE_SIGN                   ("QS-ISOTOPE-SIGN", "Isotope signal", "Isotope signal", Double.class),
    QUALSCORE_ISOTOPE_NOISE                  ("QS-ISOTOPE-NOISE", "Isotope noise", "Isotope noise", Double.class),
    QUALSCORE_ISOTOPE_STN                    ("QS-ISOTOPE-STN", "Isotope S/N ratio", "Isotope S/N ratio", Double.class),

    QUALSCORE_NUM_PEAKS                      ("QS-NUM-PEAKS","Number of peaks","Number of peaks",Integer.class),
    QUALSCORE_AVG_BY_INTENSITY               ("QS-AVG_INTENSITY", "Average intensity per peak", "Average intensity per peak", Double.class),
    QUALSCORE_STD_INTENSITY                  ("QS-STD-INTENSITY", "Standard deviation of intensities", "Standard deviation of intensitie", Double.class),
    QUALSCORE_MZ_95_INTENSITY                ("QS-MZ-95-INTNESITY", "m/z range containing 95% of intensity", "m/z range containing 95% of intensity", Double.class),
    QUALSCORE_MZ_50_INTENSITY                ("QS-MZ-50-INTNESITY", "m/z range containing 50% of intensity", "m/z range containing 50% of intensity", Double.class),
    QUALSCORE_TIC_MZ                         ("QS-TIC-MZ","TIC per m/z","TIC per m/z",Double.class),
    QUALSCORE_MASS_GAP                       ("QS-MASS-GAP","SD mass-gap consec. peaks","SD mass-gap consec. peaks",Double.class),
    QUALSCORE_NEIGHGOR_2DA                   ("QS-NEIGHGOR-2DA","Average number of neighbor peaks within 2 Da","Average number of neighbor peaks within 2 Da",Double.class)
    ;
    private String code;
    private String title;
    private String fullTitle;
    private Class type;

    private SpectrumFeatureType(String code, String title, String fullTitle, Class type) {

        this.code = code;
        this.title = title;
        this.fullTitle = fullTitle;
        this.type = type;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public String getFullTitle() {
        return fullTitle;
    }

    public void setFullTitle(String fullTitle) {
        this.fullTitle = fullTitle;
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
