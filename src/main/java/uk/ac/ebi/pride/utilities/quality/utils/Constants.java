package uk.ac.ebi.pride.utilities.quality.utils;

/**
 * Created with IntelliJ IDEA.
 * @author ypriverol
 * To change this template use File | Settings | File Templates.
 */
public class Constants {

    /**
     * The masses of all amino acids (duplicates don't have to be listed multiple times).
     */
    public static final double[] AAMasses =
            new double[] {
                    57.05,  // Gly
                    71.08,  // Ala
                    87.08,  // Ser
                    97.12,  // Pro
                    99.13,  // Val
                    101.11, // Thr
                    103.14, // Cys
                    113.16, // Leu and Ile
                    114.10, // Asparagine
                    115.09, // Aspartic acid
                    128.13, // Glutamine
                    128.62, // Lysine
                    129.12, // Glutamic acid
                    131.19, // Methionine
                    137.14, // Histidine
                    147.18, // Phenylalanine
                    156.19, // Arginine
                    163.18, // Tyrosine
                    186.21, // Tryptophan
                    99.03   // Valine
            };

    public static final double WINDOW = 1.5;

    public static final int maxPeaksPer1000Da = 400;

}
