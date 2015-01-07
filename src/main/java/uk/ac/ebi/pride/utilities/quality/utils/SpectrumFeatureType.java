package uk.ac.ebi.pride.utilities.quality.utils;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;

import java.util.ArrayList;
import java.util.List;

/**
 * Spectrum Features
 */
public enum SpectrumFeatureType {

    XXArea                                   ("XXArea",           "compute the XArea of an Spectrum","compute the XArea of an Spectrum", Double.class),
    QUALSCORE_AA_MASS_DIFF                   ("QSMASS",           "AA mass differences, number of", "Plain number of peak-distances that could be amino acids (normalized by dividing this number by the square of the total number of peaks)",Double.class),
    QUALSCORE_AA_MASS_ABUNDANCE_WEIGHTED     ("QSMASSWT",         "AA mass differences, abundance-weighted",  "Feature 4 is roughly equivalent to feature 1 (number of AA distances) but it is an abundance-weighted version (FR)", Double.class),
    QUALSCORE_AA_MASS_TAG_AVG_LONGER         ("QSTAG",            "Sequence tag, average longest",      "(Average) Based on chains of consecutive amino-acid-like distances: for every peak <I>i</I> this <sub>i</sub> will represent the maximum number of steps of amino-acid distance can be taken towards lower m/z values",Double.class),
    QUALSCORE_AA_MASS_TAG_LONGER             ("QSTAGLNG",         "Sequence tag, longest", "(Longest) Based on chains of consecutive amino-acid-like distances: for every peak <I>i</I> this <sub>i</sub> will represent the maximum number of steps of amino-acid distance can be taken towards lower m/z values",Double.class),

    QUALSCORE_Complement_B_Y_IONS_3vs3       ("QSCOMPL3V3",       "Parent mass tolerance 0.3 versus 3.0", "Parent mass tolerance 0.3 versus 3.0", Double.class),
    QUALSCORE_Complement_B_Y_IONS_1vs3       ("QSCOMPL1V3",       "Parent mass tolerance 1.0 versus 3.0", "Parent mass tolerance 1.0 versus 3.0", Double.class),
    QUALSCORE_Complement_B_Y_IONS_iso1       ("QSCOMPLISO1",      "Isotopes/background signal, parent_mass+1/parent_mass-10", "Isotopes/background signal, parent_mass+1/parent_mass-10", Double.class),
    QUALSCORE_Complement_B_Y_IONS_sgn        ("QSCOMPLSIGN",      "Complements4 (main signal/background, PM+0/PM-10)", "Complements4 (main signal/background, PM+0/PM-10)", Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_SC ("QSCOMPL3CHARGE",   "Complements: simple count","Complements: simple count",Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_AB ("QSCOMPL3CHARGEAB", "Complements: divide by averaged background","Complements: divide by averaged background",Double.class),
    QUALSCORE_Complement_B_Y_IONS_3CHARGE_SA ("QSCOMPL3CHARGESA", "Complements: subtract averaged background","Complements: subtract averaged background",Double.class),

    QUALSCORE_CROSSCORR_B_Y_IONS             ("QSCROSSBY",        "An independ identification-independent cross-correlation score", "An independ identification-independent crosscorrelation score", Double.class),

    QUALSCORE_AMONIA_17_SC                   ("QSAMONIA17",       "Neutral losses: ammonia (-17) simple count","Neutral losses: ammonia (-17) simple count",Double.class),
    QUALSCORE_AMONIA_17_AVGB                 ("QSAMONIA17AV",     "Neutral losses: ammonia (-17) divided by averaged background","Neutral losses: ammonia (-17) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_17_AVGMB                ("QSAMONIA17AVG",    "Neutral losses: ammonia (-17) count minus averaged background","Neutral losses: ammonia (-17) count minus averaged background",Double.class),
    QUALSCORE_AMONIA_18_SC                   ("QSAMONIA18",       "Neutral losses: water (-18) simple count","Neutral losses: water (-18) simple count",Double.class),
    QUALSCORE_AMONIA_18_AVGB                 ("QSAMONIA18AV",     "Neutral losses: ammonia (-18) divided by averaged background","Neutral losses: ammonia (-18) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_18_AVGMB                ("QSAMONIA18AVG",    "Neutral losses: ammonia (-18) count minus averaged background","Neutral losses: ammonia (-18) count minus averaged background",Double.class),
    QUALSCORE_AMONIA_28_SC                   ("QSAMONIA28",       "Neutral losses: water (-28) simple count","Neutral losses: water (-28) simple count",Double.class),
    QUALSCORE_AMONIA_28_AVGB                 ("QSAMONIA28AV",     "Neutral losses: ammonia (-28) divided by averaged background","Neutral losses: ammonia (-28) divided by averaged background",Double.class),
    QUALSCORE_AMONIA_28_AVGMB                ("QSAMONIA28AVG",    "Neutral losses: ammonia (-28) count minus averaged background","Neutral losses: ammonia (-28) count minus averaged background",Double.class),

    QUALSCORE_ISOTOPE_SIGN                   ("QSISOTOPESIGN",    "Isotope signal", "Isotope signal", Double.class),
    QUALSCORE_ISOTOPE_NOISE                  ("QSISOTOPENOISE",   "Isotope noise", "Isotope noise", Double.class),
    QUALSCORE_ISOTOPE_STN                    ("QSISOTOPESTN",     "Isotope S/N ratio", "Isotope S/N ratio", Double.class),

    QUALSCORE_NUM_PEAKS                      ("QSPEAKSNUM",       "Number of peaks","Number of peaks",Integer.class),
    QUALSCORE_AVG_BY_INTENSITY               ("QSINTENSITYAVG",   "Average intensity per peak", "Average intensity per peak", Double.class),
    QUALSCORE_STD_INTENSITY                  ("QSSTDINTENSITY",   "Standard deviation of intensities", "Standard deviation of intensitie", Double.class),
    QUALSCORE_MZ_95_INTENSITY                ("QSINTNESITYMZ95",  "m/z range containing 95% of intensity", "m/z range containing 95% of intensity", Double.class),
    QUALSCORE_MZ_50_INTENSITY                ("QSINTNESITYMZ50",  "m/z range containing 50% of intensity", "m/z range containing 50% of intensity", Double.class),
    QUALSCORE_TIC_MZ                         ("QSTICMZ",          "TIC per m/z","TIC per m/z",Double.class),
    QUALSCORE_MASS_GAP                       ("QSMASSGAP",        "SD mass-gap consec. peaks","SD mass-gap consec. peaks",Double.class),
    QUALSCORE_NEIGHGOR_2DA                   ("QSNEIGHGOR2DA",    "Average number of neighbor peaks within 2 Da","Average number of neighbor peaks within 2 Da",Double.class);


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
