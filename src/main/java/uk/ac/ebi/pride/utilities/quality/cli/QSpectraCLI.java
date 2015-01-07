package uk.ac.ebi.pride.utilities.quality.cli;

import org.apache.commons.cli.*;
import uk.ac.ebi.pride.utilities.quality.utils.ProcessingType;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

/**
 * Define a set of properties to be run and compute them, this is important due the research an computational
 * time. It is interesting to run some particular features rather than all of them.
 *
 */
public class QSpectraCLI {

    private static Options options;

    private static void initOptions() {

        options = new Options();

        options.addOption("h", "help", false, "show help.");

        for(SpectrumFeatureType feature: SpectrumFeatureType.values())
            options.addOption(feature.getCode(),false,feature.getTitle());


        for(ProcessingType processing: ProcessingType.values())
            options.addOption(processing.getCode(), false, processing.getTitle());

        options.addOption("out", "output",true,"tab separated file output");

        options.addOption("inClusterMgf","input-cluster-mgf", true, "input file");

        options.addOption("inMzid", "in-mzid", true, "mzindetml input file");

        options.addOption("inRelatedSpectra", true, "spectra related with mzid file");

        options.addOption("allFeatures", false, "compute all features for each spectra");
    }

    public static void main(String[] args) throws ParseException {

        initOptions();

        parse(args);
    }

    public static void parse(String[] args){

        CommandLine cmd = null;
        CommandLineParser parser = new BasicParser();


        try {
            cmd = parser.parse(options, args);
            if (cmd.hasOption("h"))
                help(options);
            if(!cmd.hasOption("inMzid") && !cmd.hasOption("inClusterMgf")){
                help(options);
            }
        } catch (ParseException e) {
                help(options);
        }
    }

    /**
     * Print the Help
     * @param options
     */
    private static void help(Options options) {
        HelpFormatter formater = new HelpFormatter();
        formater.printHelp("Main", options);
        System.exit(0);
    }
}
