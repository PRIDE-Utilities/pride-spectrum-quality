package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.PeakControllerImpl;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.io.File;
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
        for(Comparable id: mzIdentMLController.getSpectrumIds()){
            Spectrum spectrum = mzIdentMLController.getSpectrumById(id);
            if(DataAccessUtilities.getMsLevel(spectrum) > 1 && spectrum.getPrecursors() != null){
                System.out.print("\nid:\t" + id);
                Collection<Map<SpectrumFeatureType, Object>> features = generator.computeFeatureForSpectrum(spectrum, DataAccessUtilities.getPrecursorCharge(spectrum.getPrecursors())).values();
                for(Map<SpectrumFeatureType, Object> featuresValues: features)
                    for(SpectrumFeatureType spectrumFeatureType: featuresValues.keySet())
                        System.out.print("\tfeature:\t" + spectrumFeatureType.getCode() + "\t:\t" + featuresValues.get(spectrumFeatureType));
            }
        }

    }
}