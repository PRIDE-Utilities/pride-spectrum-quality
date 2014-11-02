package uk.ac.ebi.pride.utilities.quality.utils.svm;

import java.util.Vector;
import weka.core.Utils;

/**
 * Class to Normalize an Array of doubles, also can be use to extract
 * the Max and Min value in a Double Array.
 * User: yperez, enriquea
 * Date: 2/08/13
 */

public class NormalizeVector {

    /*
     * This method normalize a vector based on max element
     */
    public static double[] normalizeMax (double[] dbs, double  norm){
        Utils.normalize(dbs, norm);
        return dbs;
    }

    /*
     * This method normalize a vector based on mean
     */
    public static double[] normalizeMean (double[] dbs, double norm){
        Utils.normalize(dbs, norm);
        return dbs;
    }

    /*
     * This method return the max element
     */

    public static double getMaxElement (double[] dbs){
        return dbs[Utils.maxIndex(dbs)];
    }

    /*
     * This method return the mean vector
     */
    public static double getMean (double[] dbs){
        return Utils.mean(dbs);
    }

    /*
     * Converter vector to array []
     */
    private double [] toArrayDouble (Vector <Double> vector){
        double [] dbs = new double[vector.size()];
        int i = 0;
        for(double db: vector){
            dbs[i]= db;
            i++;
        }
      return dbs;
    }
}
