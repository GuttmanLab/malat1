package edu.caltech.lncrna.malat1.smit.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.PairedEndAlignment;
import edu.caltech.lncrna.bio.alignment.SingleReadAlignment;
import edu.caltech.lncrna.bio.annotation.BedFileRecord;
import edu.caltech.lncrna.bio.annotation.Block;
import edu.caltech.lncrna.bio.annotation.Strand;
import edu.caltech.lncrna.bio.datastructures.GenomeTree;
import edu.caltech.lncrna.bio.io.BedParser;
import edu.caltech.lncrna.bio.io.PairedEndBamParser;

public class PolymeraseThreePrimeCounter {

    private Path geneBedPath;
    private Path bamPath;
    private Path outputPath;
    
    private final GenomeTree<Block> introns = new GenomeTree<>();
    private final GenomeTree<PairedEndAlignment> fragments =
            new GenomeTree<>();
    private final Map<Splicing, Map<Integer, Integer>> countsByPosition =
            new HashMap<>();
    
    private static int MAX_INSERT_SIZE = 1000;
    
    public static void main(String[] args) throws IOException {
        (new PolymeraseThreePrimeCounter())
            .parseArgs(args)
            .loadIntrons()
            .loadReads()
            .calculatePolPositionsForIntrons()
            .printCounts();
    }
  
    /**
     * Parses the command-line arguments.
     * <p>
     * Any parse errors cause the program to print a help menu, then exit with
     * a -1 status code.
     * @param args - the command-line arguments
     */
    private PolymeraseThreePrimeCounter parseArgs(String[] args) {
        Options options = new Options();

        options.addOption(Option.builder()
                .longOpt("bam")
                .desc("input BAM file")
                .argName("BAM")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("genes")
                .desc("BED file of genes or transcripts")
                .argName("BED")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("output")
                .desc("output file")
                .argName("OUT")
                .hasArg()
                .required()
                .build());

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);
            bamPath = Paths.get(cmd.getOptionValue("bam"));
            geneBedPath = Paths.get(cmd.getOptionValue("genes"));
            outputPath = Paths.get(cmd.getOptionValue("output"));
        } catch (ParseException e) {
            printHelp(options);
            System.exit(-1);
        }
        
        return this;
    }
    
    /**
     * Creates and prints a help menu to the console.
     * @param opts - the options to include in the help menu
     */
    private void printHelp(Options opts) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setDescPadding(0);
        String header = "\n";
        String footer = "\n";
        formatter.printHelp("java -jar PolymeraseThreePrimeCounter.jar", header,
                opts, footer, true);
    }
    
    /**
     * Loads records from a BED file into memory. Only keeps those records
     * which do not overlap any other record.
     */
    private PolymeraseThreePrimeCounter loadIntrons() {
        GenomeTree<BedFileRecord> genes = new GenomeTree<>();
        
        try (BedParser p = new BedParser(geneBedPath)) {
            p.stream().forEach(x -> genes.insert(x));
        }
        
        try (BedParser p = new BedParser(geneBedPath)) {
            p.stream()
             .filter(x -> genes.getNumOverlappers(x) <= 1)
             .flatMap(x -> x.getIntronBlockStream())
             .forEach(x -> introns.insert(x));
        }
        
        return this;
    }
    
    /**
     * Loads fragments from a paired-end BAM file into memory. Only keeps those
     * fragments which (hull-)overlap a gene in <code>genes</code>. 
     */
    private PolymeraseThreePrimeCounter loadReads() {

        try (PairedEndBamParser records = new PairedEndBamParser(bamPath)) {
            records.getAlignmentStream()
                   .filter(x -> introns.overlaps(x.getHull()))
                   .forEach(x -> fragments.insert(x));
        }

        return this;
    }
    
    private PolymeraseThreePrimeCounter calculatePolPositionsForIntrons() {
        
        countsByPosition.put(Splicing.SPLICED, new TreeMap<>());
        countsByPosition.put(Splicing.UNSPLICED, new TreeMap<>());
        
        for (Block intron : introns) {
            calculatePolPositionsForIntron(intron);
        }
        return this;
    }
    
    private void calculatePolPositionsForIntron(Block intron) {
        Iterator<PairedEndAlignment> overlappingFragments = 
                fragments.getHullOverlappers(intron);
        
        while (overlappingFragments.hasNext()) {
            PairedEndAlignment fragment = overlappingFragments.next();
            calculatePolPositionForFragmentAndIntron(fragment, intron);
        }
    }
    
    /**
     * Calculates the Pol II position implied by a paired-end fragment relative
     * to a given intron and increments the appropriate count.
     * @param fragment - the paired-end alignment
     * @param intron - the intron
     */
    private void calculatePolPositionForFragmentAndIntron(
            PairedEndAlignment fragment, Block intron) {

        Splicing splice = getSplicingRelationship(fragment, intron);
        
        if (splice.isValid) {
            int intronThreePrime = intron.getThreePrimePosition();
            int polymerasePosition = fragment.getThreePrimePosition();
            int relativePolPosition = fragment.getStrand().equals(Strand.POSITIVE)
                    ? polymerasePosition - intronThreePrime
                    : intronThreePrime - polymerasePosition;
            Map<Integer, Integer> spliceMap = countsByPosition.get(splice);
            int count = spliceMap.getOrDefault(relativePolPosition, 0);
            spliceMap.put(relativePolPosition, count + 1);
        }
    }
    
    /**
     * Gets the splicing state of the region that a given read-pair aligns to.
     * <p>
     * Determines if the read-pair is spliced, unspliced, etc. with respect to
     * the given intron. Ignores other introns that the read-pair may overlap.
     * @param fragment - the paired-end fragment
     * @param intron - the intron
     */
    private Splicing getSplicingRelationship(PairedEndAlignment fragment,
            Block intron) {

        if (fragment.getInsertSize() > MAX_INSERT_SIZE) {
            return Splicing.SPLICED;
        }
        
        SingleReadAlignment read1 = fragment.getFirstReadInPair();
        Region region1 = getAlignmentRegion(read1, intron);
        
        SingleReadAlignment read2 = fragment.getSecondReadInPair();
        Region region2 = getAlignmentRegion(read2, intron);
        
        return region1.with(region2);
    }
    
    /**
     * Gets the type of region that a given read aligns to, relative to a given
     * intron
     * @param read - the aligned read
     * @param intron - the intron
     */
    private Region getAlignmentRegion(SingleReadAlignment read, Block intron) {

        // Read is entirely within the intron
        if (intron.contains(read)) {
            return Region.INTRON;
        }
        
        // Read not entirely within the intron, but overlaps it
        // The read must span splice exon-intron splice junction
        if (intron.overlaps(read)) {
            return Region.EXON_INTRON_JUNCTION;
        }

        // Read spans the intron with no overlap. The intron must have been
        // spliced out.
        if (read.getHull().contains(intron)) {
            return Region.EXON_EXON_JUNCTION;
        }
        
        // All other cases indicate the read is exonic.
        return Region.EXON;
    }
    
    /**
     * Prints the counts data to the output file.
     * <p>
     * File format is tab-delimited with three columns:
     * <li>PolII position relative to intron three-prime end
     * <li>Number of fragments corresponding to the position in column 1.
     * <li>The splicing state ("SPLICED" or "UNSPLICED")
     * @throws IOException
     */
    private void printCounts() throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
            printCountsForSplicing(Splicing.SPLICED, writer);
            printCountsForSplicing(Splicing.UNSPLICED, writer);
        }
    }
    
    /**
     * Prints the counts data for the given <code>Splicing</code>
     * @param splicing - the splice-state
     * @param out - the output writer
     * @throws IOException
     */
    private void printCountsForSplicing(Splicing splicing, BufferedWriter out)
            throws IOException {

        Map<Integer, Integer> spliceCounts =
                countsByPosition.getOrDefault(splicing, Collections.emptyMap());
        for (Map.Entry<Integer, Integer> entry : spliceCounts.entrySet()) {
            out.write(String.valueOf(entry.getKey()));
            out.write("\t");
            out.write(String.valueOf(entry.getValue()));
            out.write("\t");
            out.write(splicing.toString());
            out.newLine();
        }
    }

    /**
     * An enumeration of regions where a read may align relative to an intron.
     */
    private enum Region {
        
        EXON(Splicing.UNKNOWN),
        INTRON(Splicing.UNSPLICED),
        EXON_EXON_JUNCTION(Splicing.SPLICED),
        EXON_INTRON_JUNCTION(Splicing.UNSPLICED);
        
        private Splicing splice;
        
        private Region(Splicing splice) {
            this.splice = splice;
        }
        
        /**
         * Gets the splicing state implied by this region.
         * @param other - the other region
         */
        public Splicing with(Region other) {
            return splice.with(other.splice);
        }
    }
    
    /**
     * An enumeration of splice-states.
     * <p>
     * A read can map to a region that is definitely spliced, a region that is
     * definitely unspliced, or a region about which we cannot discern the
     * splice state. An invalid splicing indicates that we get conflicting
     * splicing information from the two reads in a read-pair. 
     */
    private enum Splicing {

        SPLICED(true) {
            public Splicing with(Splicing other) {
                switch (other) {
                case SPLICED:
                    return Splicing.SPLICED;
                case UNSPLICED:
                    return Splicing.INVALID;
                case UNKNOWN:
                    return Splicing.SPLICED;
                case INVALID:
                    return Splicing.INVALID;
                default:
                    throw new IllegalStateException("Default block should " +
                            "not be reached.");
                }
            }
        },
        
        UNSPLICED(true) {
            public Splicing with(Splicing other) {
                switch (other) {
                case SPLICED:
                    return Splicing.INVALID;
                case UNSPLICED:
                    return Splicing.UNSPLICED;
                case UNKNOWN:
                    return Splicing.UNSPLICED;
                case INVALID:
                    return Splicing.INVALID;
                default:
                    throw new IllegalStateException("Default block should " +
                            "not be reached.");   
                }
            }
            
        },
        
        UNKNOWN(false) {
            public Splicing with(Splicing other) {
                return other;
            }
        },
        
        INVALID(false) {
            public Splicing with(Splicing other) {
                return Splicing.INVALID;
            }
        };
        
        private boolean isValid;
        
        private Splicing(boolean isValid) {
            this.isValid = isValid;
        }
        
        /**
         * Returns the splice-state implied by this and another
         * <code>Splicing</code>.
         * <p>
         * Each read in a paired-end fragment contains (potentially
         * different) information about whether it is spliced. This method
         * returns the splice-state of the fragment implied by both of its
         * reads.
         * @param the other splicing state
         */
        public abstract Splicing with(Splicing other);
    }
}