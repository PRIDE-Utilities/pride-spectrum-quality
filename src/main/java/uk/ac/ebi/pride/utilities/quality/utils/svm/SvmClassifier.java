package uk.ac.ebi.pride.utilities.quality.utils.svm;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMOreg;
import weka.classifiers.functions.supportVector.Kernel;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SerializationHelper;


public class SvmClassifier implements Serializable{

    private SMOreg svm;
    /*An SVM instance*/

    private Evaluation eval;

    private Instances trainData;

    private double rmsd;

    //constructor
    public SvmClassifier(Instances trainData, String[] options, Kernel kernel) throws Exception{
        this.trainData = trainData;
        this.svm = new SMOreg();
        this.svm.setOptions(options);       //set the options
        this.svm.setKernel(kernel);
        this.svm.buildClassifier(trainData);    // build classifier
    }

    /*
     * This method return standard deviation obtained in the stage test the classifier
     * it is used like selection criteria.
     */
    //
    public double evaluateRMSE(Instances test) throws Exception{
        this.eval = new Evaluation(trainData);
        eval.evaluateModel(svm, test);
        double rmsd = eval.errorRate();
        return rmsd;
    }

    public double evaluateCorrelation(Instances test) throws Exception{
        this.eval = new Evaluation(trainData);
        eval.evaluateModel(svm, test);
        double rsquared = eval.correlationCoefficient();
        return rsquared;
    }


    public Instances classifingInstance(Instances dataset) throws Exception{
        dataset.setClassIndex(-1); //data set without class attribute
        dataset.setClassIndex(dataset.numAttributes() - 1);
        Instances labeled = new Instances(dataset);
        for (int i = 0; i < dataset.numInstances(); i++) {
            double clsLabel = svm.classifyInstance(dataset.instance(i));
            labeled.instance(i).setClassValue(clsLabel);
        }
        return labeled;
    }

    public double[] getCoefficients(){
        //Todo: Basically this function must be implemented using a Weka function but do not exist
        String svmClassifier = svm.toString();
        List<Double> coefficents = new ArrayList<Double>();
        double[] coeficcientsArray = null;
        if(svmClassifier.contains("weights (not support vectors):")){
          String[] svmClassStrings = svmClassifier.split("\t");
          int i = 0;
          boolean lastFound = false;
          while( i < svmClassStrings.length && !lastFound){
              String currentString = svmClassStrings[i];
              if(currentString.contains("weights (not support vectors):")){
                  i++;
                  while(i < svmClassStrings.length && !lastFound){
                      currentString = svmClassStrings[i];
                      String[] parts = currentString.split("\\s+");
                      if(parts.length < 5){
                          lastFound = true;
                      }else{
                          String number = parts[1] + parts[2];
                          coefficents.add(Double.parseDouble(number));
                      }
                      i++;
                  }
              }
              i++;
          }
        }

        if(!coefficents.isEmpty()){
            coeficcientsArray = new double[coefficents.size()];
            for(int i =0 ; i < coefficents.size(); i++) coeficcientsArray[i] = coefficents.get(i);
        }
        return coeficcientsArray;
    }


    public double classifyInstance(Instance inst) throws Exception{
        return this.eval.evaluateModelOnce(svm,inst);
    }

    public static void saveClassifier(String fileSaveName, SvmClassifier svmClassifier) throws Exception {
        SerializationHelper.write(fileSaveName, svmClassifier);
    }

    public static SvmClassifier loadClassifier(String idModel) throws Exception {
        return (SvmClassifier) SerializationHelper.read(idModel);
    }

    @Override
    public String toString() {
        return "SvmClassifier{" +
                "svm=" + svm.toString() +
                ", eval=" + eval.toString() +
                ", trainData=" + trainData.toString() +
                ", rmsd=" + rmsd +
                '}';
    }
}
