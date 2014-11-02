package uk.ac.ebi.pride.utilities.quality.utils.svm;

import weka.classifiers.functions.supportVector.PolyKernel;
import weka.classifiers.functions.supportVector.RBFKernel;
import weka.core.Instances;


public class SvmKernel {

    private static volatile RBFKernel rbfKernel = null;
    /* Radial Kernel */

    private static volatile PolyKernel polyKernel = null;
    /* Polynomial Kernel */

    private static volatile PolyKernel linealKernel = null;
    /* Lineal Kernel is an special case of Polynomial Kernel*/

    public static RBFKernel getInstanceRBFKernel(Instances trainData, int cache, double gamma) throws Exception {
        if (rbfKernel == null) {
            synchronized (RBFKernel .class){
                if (rbfKernel == null) {
                    rbfKernel = new RBFKernel(trainData,cache,gamma);
                }
            }
        }
        return rbfKernel;
    }

    public static PolyKernel getInstancePolyKernel(Instances trainData, int cache, double exp, boolean order) throws Exception {
        if (polyKernel == null) {
            synchronized (PolyKernel .class){
                if (polyKernel == null) {
                    polyKernel = new PolyKernel(trainData,cache,exp,order);
                }
            }
        }
        return polyKernel;
    }

    public static PolyKernel getInstanceLinearKernel(Instances trainData, int cache) throws Exception {
        if (linealKernel == null) {
            synchronized (PolyKernel .class){
                if (linealKernel == null) {
                    linealKernel = new PolyKernel(trainData,cache,1.0,false);
                }
            }
        }
        return linealKernel;
    }

}
