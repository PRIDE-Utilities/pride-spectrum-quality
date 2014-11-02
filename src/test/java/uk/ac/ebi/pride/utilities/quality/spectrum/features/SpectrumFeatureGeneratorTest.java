package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.PeakControllerImpl;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.io.File;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class SpectrumFeatureGeneratorTest {

    SpectrumFeatureGenerator generator;

    PeakControllerImpl mgfController;

    MzIdentMLControllerImpl mzIdentMLController;

    PrintStream outFile = null;



    @org.junit.Before
    public void setUp() throws Exception {
        URL url = SpectrumFeatureGeneratorTest.class.getClassLoader().getResource("small.mzid");
        if (url == null) {
            throw new IllegalStateException("no file for input found!");
        }
        File inputFile = new File(url.toURI());
        mzIdentMLController = new MzIdentMLControllerImpl(inputFile, true);

        URL urlmgf = SpectrumFeatureGeneratorTest.class.getClassLoader().getResource("small.mgf");
        File filems = new File(urlmgf != null ? urlmgf.getFile() : null);
        List<File> fileList = new ArrayList<File>();
        fileList.add(filems);
        mzIdentMLController.addMSController(fileList);
        outFile = new PrintStream(new File("variables.txt"));
    }

    @org.junit.After
    public void tearDown() throws Exception {
      mzIdentMLController.close();
    }

    @org.junit.Test
    public void testGetInstance() throws Exception {

    }

    @org.junit.Test
    public void testGetFeatureCount() throws Exception {

    }

    @org.junit.Test
    public void testComputeFeatureForSpectrum() throws Exception {
        generator = SpectrumFeatureGenerator.getInstance();
        int count = 0;
        for(Comparable id: mzIdentMLController.getSpectrumIds()){
            Spectrum spectrum = mzIdentMLController.getSpectrumById(id);
            if(DataAccessUtilities.getMsLevel(spectrum) > 1 && spectrum.getPrecursors() != null){
                Map<Integer, Map<SpectrumFeatureType, Object>> features = generator.computeFeatureForSpectrum(spectrum, DataAccessUtilities.getPrecursorCharge(spectrum.getPrecursors()));
                if(count == 0){
                    outFile.print("id\tidentified\t");
                    for(Integer featuresValues: features.keySet())
                        for(SpectrumFeatureType spectrumFeatureType: features.get(featuresValues).keySet())
                            outFile.print(spectrumFeatureType + "["+ featuresValues+ "]" +"\t");
                }
                outFile.println();
                outFile.print(id + "\t" + mzIdentMLController.isIdentifiedSpectrum(id) + "\t");
                Collection<Map<SpectrumFeatureType, Object>> featuresCollection = features.values();
                for(Map<SpectrumFeatureType, Object> featuresValues: featuresCollection)
                    for(SpectrumFeatureType spectrumFeatureType: featuresValues.keySet())
                        outFile.print(featuresValues.get(spectrumFeatureType) + "\t");
            }
            count++;
        }

    }
}