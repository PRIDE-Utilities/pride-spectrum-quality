package uk.ac.ebi.pride.utilities.quality.spectrum.common;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;

import java.util.*;

/**
 * Spectrum Utilities  contains a set of functions to sort peaks by Mass, by Peak intensity
 * Also it contains functions to select peaks, etc.
 * @author ypriverol
 */
public class SpectrumUtils {

    /**
     * Sort the peaks of mass spectrum by Mass. The revertOrder guaranty of the function is
     * computed from largest (revertOrder = true) to smallest or smallest to largest ((revertOrder = false)).
     * @param spectrum Spectrum Object
     * @param revertOrder
     * @return
     */
    public static Spectrum sortSpectrumByMass(Spectrum spectrum, boolean revertOrder){

        double[] intensities = spectrum.getIntensityBinaryDataArray().getDoubleArray();
        double[] masses      = spectrum.getMzBinaryDataArray().getDoubleArray();
        for(int i = 0; i < masses.length; i++){
            for(int j = i + 1; j < masses.length; j++){
                if(!revertOrder){
                    if(masses[i] > masses[j]){
                        double mass = masses[i];
                        masses[i]   = masses[j];
                        masses[j]   = mass;
                        double intensity = intensities[i];
                        intensities[i]   = intensities[j];
                        intensities[j]   = intensity;
                    }
                }else{
                    if(masses[i] < masses[j]){
                        double mass = masses[i];
                        masses[i]   = masses[j];
                        masses[j]   = mass;
                        double intensity = intensities[i];
                        intensities[i]   = intensities[j];
                        intensities[j]   = intensity;
                    }
                }
            }
        }
        spectrum.getMzBinaryDataArray().setDoubleArray(masses);
        spectrum.getIntensityBinaryDataArray().setDoubleArray(intensities);
        return spectrum;

    }

    /**
     * Sort the peaks of mass spectrum by Intensity. The revertOrder guaranty of the function is
     * computed from largest (revertOrder = true) to smallest or smallest to largest ((revertOrder = false)).
     * @param spectrum Spectrum Object
     * @param revertOrder
     * @return
     */
    public static Spectrum sortSpectrumByIntensity(Spectrum spectrum, boolean revertOrder){

        double[] intensities = spectrum.getIntensityBinaryDataArray().getDoubleArray();
        double[] masses      = spectrum.getMzBinaryDataArray().getDoubleArray();
        for(int i = 0; i < intensities.length; i++){
            for(int j = i + 1; j < intensities.length; j++){
                if(!revertOrder){
                    if(intensities[i] > intensities[j]){
                        double mass = masses[i];
                        masses[i]   = masses[j];
                        masses[j]   = mass;
                        double intensity = intensities[i];
                        intensities[i]   = intensities[j];
                        intensities[j]   = intensity;
                    }
                }else{
                    if(intensities[i] < intensities[j]){
                        double mass = masses[i];
                        masses[i]   = masses[j];
                        masses[j]   = mass;
                        double intensity = intensities[i];
                        intensities[i]   = intensities[j];
                        intensities[j]   = intensity;
                    }
                }
            }
        }
        spectrum.getMzBinaryDataArray().setDoubleArray(masses);
        spectrum.getIntensityBinaryDataArray().setDoubleArray(intensities);
        return spectrum;

    }

    public static double getTotIonCurrentCount(Spectrum spectrum){
        double ionCurrentCount = 0.0;
        for(int i = 0; i < spectrum.getIntensityBinaryDataArray().getDoubleArray().length; i++)
            ionCurrentCount = ionCurrentCount + spectrum.getIntensityBinaryDataArray().getDoubleArray()[i];
        return ionCurrentCount;
    }




}
