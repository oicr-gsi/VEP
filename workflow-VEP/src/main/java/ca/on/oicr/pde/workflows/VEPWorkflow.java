package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.HashMap;
import java.util.Map;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class VEPWorkflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;

    // Input Data
    private String inputVCF;
    private String inputVCFindex;
    private String outputFilenamePrefix;

    //vcf2maf
    private String vcf2mafpl;
    private String perlPath;
    private String vcf2mafPath;
    private String perl;

    // preprocess
    String bgzip;
    String tabix;
    String tabixVersion;
    String bedtools;
    String bcftools;
    String vcftools;

    //params
    private String hgBuild;
    private String species;
    private double vafFilter = 0.7;
    private Integer bufferSize = 200;
    private Integer acFilter = 10;
    private String vepPath;
    private String vepData;
    private String retainInfo;
    private String oncoKBPath;

    // additional params
    private String additionalParams;

    // path var
    String pathfixExport;

    //Memory allocation
    private Integer VEPMem;

    //ref Data
    private String refFasta;
    private String exacVCF;

    // target bed file
    private String targetBedFile;
    // tglFrequency file
    private String freqTextFile;

    private boolean manualOutput;
    private String queue;

//    // metatypes
    private final String TXT_GZ_METATYPE = "application/txt-gz";
    private final String VCF_GZ_METATYPE1 = "application/vcf-4-gzip";
    private final String VCF_GZ_METATYPE2 = "application/vcf-gz";
    private final String VCF_METATYPE = "application/vcf";
    private final String VCF_TBI_METATYPE = "application/tbi";

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            inputVCF = getProperty("input_vcf_file");
            inputVCFindex = inputVCF + ".tbi";
            outputFilenamePrefix = getProperty("output_filename_prefix");
//            normalSamplePrefix = getOptionalProperty("matched_normal_name", "matched");

            // vcf2maf
            vcf2mafpl = getProperty("vcf2maf");
            perlPath = getProperty("perl_path");
            vcf2mafPath = getProperty("vcf2maf_path");
            perl = getProperty("tgl_perl");
            oncoKBPath = getProperty("oncokb");

            // preprocess
            tabixVersion = getProperty("tabix_version");
            bgzip = getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip";
            tabix = getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix";
            bcftools = getProperty("bcftools");
            bedtools = getProperty("bedtools");
            vcftools = getProperty("vcftools");

            // ref fasta
            refFasta = getProperty("ref_fasta");
            // exac vcf
            exacVCF = getProperty("exac_vcf");
            // target bed file
            targetBedFile = getProperty("target_bed");
            //tglfreq file
            freqTextFile = getOptionalProperty("freq_file", null);

            // VEP
            vafFilter = Double.parseDouble(getOptionalProperty("vaf_filter", "0.7"));
            species = getOptionalProperty("species", "homo_sapiens");
            hgBuild = getOptionalProperty("hg_version", "GRCh37");
            vepPath = getProperty("vep_path");
            vepData = getProperty("vep_data");
            bufferSize = Integer.parseInt(getOptionalProperty("buffer_size", "200"));
            acFilter = Integer.parseInt(getOptionalProperty("max_ac_filter", "10"));

            //additional params
            additionalParams = getOptionalProperty("additional_arg", null);

            //
            pathfixExport = "# perl/5.22.2-tgl\n";
            pathfixExport = pathfixExport + "export LD_LIBRARY_PATH=" + this.perlPath + "/lib:$LD_LIBRARY_PATH;\n";
            pathfixExport = pathfixExport + "export PERL5LIB=" + this.perlPath + "/lib:$PERL5LIB\n";
            pathfixExport = pathfixExport + "export PATH=" + this.perlPath + "/bin:$PATH\n";
            pathfixExport = pathfixExport + "# vep/92\n";
            pathfixExport = pathfixExport + "export PATH=" + this.vepPath + ":$PATH\n";
            pathfixExport = pathfixExport + "export PATH=" + this.vepPath + "/htslib:$PATH\n";
            pathfixExport = pathfixExport + "export PATH=" + this.vepPath + "/samtools/bin:$PATH\n";
            pathfixExport = pathfixExport + "export PERL5LIB=" + this.vepPath + ":$PATH\n";
            pathfixExport = pathfixExport + "export VEP_PATH=" + this.vepPath + ";\n";
            pathfixExport = pathfixExport + "export VEP_DATA=" + this.vepData + ";\n";
            pathfixExport = pathfixExport + "# vcf2maf\n";
            pathfixExport = pathfixExport + "export PATH=" + this.vcf2mafPath + ":$PATH;\n";
            pathfixExport = pathfixExport + "export PATH=" + this.oncoKBPath + ":$PATH;\n";
            pathfixExport = pathfixExport + "export PATH=" + this.vcftools + ":$PATH;\n";
            pathfixExport = pathfixExport + "export PATH=" + this.bedtools + ":$PATH;\n";
            pathfixExport = pathfixExport + "export PATH=" + this.bcftools + ":$PATH;\n";

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            VEPMem = Integer.parseInt(getProperty("vep_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("inVCF");
        if (inputVCF.endsWith("gz")) {
            if (inputVCF.contains("mutect2")) {
                file0.setSourcePath(inputVCF);
                file0.setType(VCF_GZ_METATYPE2);
                file0.setIsInput(true);
            } else {
                file0.setSourcePath(inputVCF);
                file0.setType(VCF_GZ_METATYPE1);
                file0.setIsInput(true);     
            }
            SqwFile file1 = this.createFile("inVCFTBI");
            file1.setSourcePath(inputVCFindex);
            file1.setType(VCF_TBI_METATYPE);
            file1.setIsInput(true);
        } else {
            file0.setSourcePath(inputVCF);
            file0.setType(VCF_METATYPE);
        }
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        Job parentJob = null;
        String inVCF = getFiles().get("inVCF").getProvisionedPath();
        String VCFtbi = getFiles().get("inVCFTBI").getProvisionedPath();
        String mafFile = this.dataDir + this.outputFilenamePrefix + ".maf.txt";

        // extract sample names first
        Job extractSampleNames = this.extractSampleNames(inVCF);
        parentJob = extractSampleNames;

        // unmatched VCFs are labelled as "tumor_only" VCFs
        if (inVCF.contains("tumor_only")) {
            Job preprocessUnmatchedVCF = handleUnmatchedVCF(inVCF);
            preprocessUnmatchedVCF.addParent(parentJob);
            parentJob = preprocessUnmatchedVCF;
            inVCF = this.tmpDir + this.outputFilenamePrefix + ".unmatched.vcf.gz";
        }

        // subset VCF 
        HashMap<String, Job> getTargetVCF = preProcessVCF(inVCF);
        String subsetVCF = getTargetVCF.keySet().toArray()[0].toString();
        Job preprocess = getTargetVCF.get(subsetVCF);
        preprocess.addParent(parentJob);
        parentJob = preprocess;

        // provision out subset VCF
        SqwFile targetVCF = createOutputFile(subsetVCF, VCF_GZ_METATYPE1, this.manualOutput);
        targetVCF.getAnnotations().put("Target_VCF", "VEP");
        parentJob.addFile(targetVCF);

        // provision out subset VCF TBI
        SqwFile targetVCFtbi = createOutputFile(subsetVCF + ".tbi", VCF_TBI_METATYPE, this.manualOutput);
        targetVCFtbi.getAnnotations().put("Target_VCF_tbi", "VEP");
        parentJob.addFile(targetVCFtbi);

        // annotate frequency
        String intVCF = subsetVCF;
        String tglFreqVCF;
        if (this.freqTextFile != null) {
            Job tglFreq = TGLFreqAnnotation(intVCF);
            tglFreqVCF = this.tmpDir + this.outputFilenamePrefix + "_final.tglfreq.vcf";
            this.retainInfo = "TGL_Freq";
            tglFreq.addParent(parentJob);
            parentJob = tglFreq;
        } else {
            tglFreqVCF = inVCF;
        }

        // run vcf to maf
        Job vcf2MAF = runVcf2Maf(tglFreqVCF, mafFile);
        vcf2MAF.addParent(parentJob);
        parentJob = vcf2MAF;

        // oncokb annotator
        Job oncoKBAnnotate = getWorkflow().createBashJob("oncokb_annotate");
        //cmd.addArgument(pathfixExport);
        oncoKBAnnotate.setCommand(pathfixExport + this.oncoKBPath
                + "/" + "MafAnnotator.py -i "
                + mafFile + " -o " + mafFile.replace(".txt", ".oncoKB.txt"));
        oncoKBAnnotate.addParent(parentJob);
        parentJob = oncoKBAnnotate;
        parentJob.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        parentJob.setQueue(getOptionalProperty("queue", ""));

        // zip maf file
        String oncoKBMafFile = mafFile.replace(".txt", ".oncoKB.txt");
        Job zipMafFile = getWorkflow().createBashJob("zip_maf");
        zipMafFile.setCommand(bgzip + " -c "
                + oncoKBMafFile + " > " + mafFile + ".gz");
        zipMafFile.addParent(parentJob);
        parentJob = zipMafFile;
        parentJob.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        parentJob.setQueue(getOptionalProperty("queue", ""));

        // Provision out maf.txt file
        SqwFile outMaf = createOutputFile(mafFile + ".gz",
                TXT_GZ_METATYPE, this.manualOutput);
        outMaf.getAnnotations().put("MAF", "VEP");
        parentJob.addFile(outMaf);

    }

    private Job runVcf2Maf(String inVCF, String outputMAF) {
        Job runVCF2MAF = getWorkflow().createBashJob("vcf2maf");
        Command cmd = runVCF2MAF.getCommand();
        cmd.addArgument(pathfixExport);
        // command to parse sample_names file
        cmd.addArgument("if [[ `cat " + this.tmpDir + "sample_names "
                + "| tr \",\" \"\\n\" | wc -l` == 2 ]]; then \n"
                + "for item in `cat " + this.tmpDir + "sample_names"
                + " | tr \",\" \"\\n\"`; do "
                + "if [[ $item == \"NORMAL\" || $item == *_R_* || $item == *BC*"
                + "  || $item == \"unmatched\" ]]; then"
                + " NORM=$item; else TUMR=$item; fi; done \n"
                + "else "
                + "TUMR=`cat " + this.tmpDir + "sample_names "
                + "| tr -d \",\"`; NORM=\"unmatched\"; fi\n\n"); //
        cmd.addArgument(this.perl + " " + this.vcf2mafpl);
        cmd.addArgument("--species " + this.species);
        cmd.addArgument("--ncbi-build " + this.hgBuild);
        cmd.addArgument("--input-vcf " + inVCF);
        cmd.addArgument("--output-maf " + outputMAF);
        cmd.addArgument("--tumor-id $TUMR");
        cmd.addArgument("--normal-id $NORM");
        cmd.addArgument("--vcf-tumor-id $TUMR");
        cmd.addArgument("--vcf-normal-id $NORM");
        cmd.addArgument("--ref-fasta " + this.refFasta);
        cmd.addArgument("--filter-vcf " + this.exacVCF);
        cmd.addArgument("--max-filter-ac " + this.acFilter);
        cmd.addArgument("--vep-path " + this.vepPath);
        cmd.addArgument("--vep-data " + this.vepData);
        if (this.freqTextFile != null) {
            cmd.addArgument("--retain-info " + this.retainInfo);
        }
        cmd.addArgument("--min-hom-vaf " + Double.toString(this.vafFilter));
        cmd.addArgument("--buffer-size " + Integer.toString(this.bufferSize));
        if (this.additionalParams != null) {
            cmd.addArgument(this.additionalParams);
        }
        runVCF2MAF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        runVCF2MAF.setQueue(getOptionalProperty("queue", ""));
        return runVCF2MAF;
    }

    // merge VCFs for mutect2 unmatched
    private Job handleUnmatchedVCF(String inVCF) {
        String tempTumorVCF = this.tmpDir + this.outputFilenamePrefix + ".vcf";
        String tempMutect2VCF = this.tmpDir + this.outputFilenamePrefix
                + ".unmatched.vcf";
        Job mergeMutect2VCF = getWorkflow().createBashJob("preprocess_unmatched");
        String inputUnmatchedVCF = this.tmpDir + this.outputFilenamePrefix
                + "_input" + ".vcf.gz";
        Command cmd = mergeMutect2VCF.getCommand();
        cmd.addArgument(pathfixExport);
        cmd.addArgument("zcat " + inVCF
                + " | sed \"s/QSS\\,Number\\=A/QSS\\,Number\\=\\./\" | "
                + this.bgzip + " -c > " + inputUnmatchedVCF + ";\n"); // fix QSS header
        cmd.addArgument(tabix + " -p vcf " + inputUnmatchedVCF + ";\n"); //tabix index the vcf 
        cmd.addArgument("if [[ `cat " + this.tmpDir + "sample_names | "
                + "tr \",\" \"\\n\" | wc -l` == 2 ]]; then \n"
                + "for item in `cat " + this.tmpDir + "sample_names"
                + " | tr \",\" \"\\n\"`; do "
                + "if [[ $item == \"NORMAL\" || $item == *_R_* ]]; "
                + "then NORM=$item; else TUMR=$item; fi; done \n"
                + "else TUMR=`cat " + this.tmpDir + "sample_names "
                + "| tr -d \",\"`; NORM=\"unmatched\"; fi\n\n");
        cmd.addArgument("echo -e \"$TUMR\\n$NORM\" > "
                + this.tmpDir + this.outputFilenamePrefix + "_header \n\n"); // create header file
        cmd.addArgument(bcftools + "/bcftools  merge " + inputUnmatchedVCF + " "
                + inputUnmatchedVCF + " --force-samples >"
                + tempTumorVCF + ";\n"); // merge VCFs
        cmd.addArgument(bcftools + "/bcftools reheader -s "
                + this.tmpDir + this.outputFilenamePrefix + "_header "
                + tempTumorVCF + ">" + tempMutect2VCF + ";\n"); //reheader merged VCF
        cmd.addArgument(bgzip + " -c " + tempMutect2VCF + " > "
                + tempMutect2VCF + ".gz" + ";\n"); // bgzip
        cmd.addArgument(tabix + " -p vcf " + tempMutect2VCF + ".gz"); // tabix index
        mergeMutect2VCF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        mergeMutect2VCF.setQueue(getOptionalProperty("queue", ""));
        return mergeMutect2VCF;
    }

    private Job TGLFreqAnnotation(String inVCF) {
        String intermediateVCF = inVCF.replace(".vcf.gz", "_temp.vcf");
        String freqAnnotVCF = this.tmpDir + this.outputFilenamePrefix
                + "_tmp_tglfreq.vcf";
        String finalTGLFreqAnnotVCF = this.tmpDir + this.outputFilenamePrefix
                + "_final.tglfreq.vcf";
        Job annotateTGLFreq = getWorkflow().createBashJob("tgl_freq");
        Command cmd = annotateTGLFreq.getCommand();
        cmd.addArgument(pathfixExport);
        cmd.addArgument(bcftools + "/bcftools annotate -a " + this.freqTextFile);
        cmd.addArgument("-c CHROM,POS,REF,ALT,TGL_Freq");
        cmd.addArgument("-h <(echo '##INFO=<ID=TGL_Freq,Number=.,"
                + "Type=Float,Description=\"Variant Frequency Among "
                + "TGL Tumours (MuTect2 Artifact Detection)\">')");
        cmd.addArgument(inVCF + " > " + intermediateVCF + ";\n");
        cmd.addArgument(this.bgzip + " -c " + intermediateVCF + " > "
                + intermediateVCF + ".gz" + ";\n");
        cmd.addArgument(this.tabix + " -p vcf " + intermediateVCF + ".gz" + ";\n");
        cmd.addArgument("echo \"Marking novel variants as TGL_Freq=0.0\"\n");
        cmd.addArgument(bcftools + "/bcftools annotate -a " + this.freqTextFile);
        cmd.addArgument("-c CHROM,POS,REF,ALT,TGL_Freq");
        cmd.addArgument("-m \"-TGL_Freq=0.0\" ");
        cmd.addArgument(intermediateVCF + ".gz" + " > " + freqAnnotVCF + ";\n");
        cmd.addArgument("cat " + freqAnnotVCF
                + " | grep -v \"Sites not listed in OICR_Freq=0.0\" > "
                + finalTGLFreqAnnotVCF + ";\n");
        annotateTGLFreq.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        annotateTGLFreq.setQueue(getOptionalProperty("queue", ""));
        return annotateTGLFreq;
    }

    private HashMap<String, Job> preProcessVCF(String inVCF) {
        HashMap<String, Job> hmap = new HashMap<String, Job>();
        int index = inVCF.lastIndexOf("/");
        String vcfName = inVCF.substring(index + 1);
        String newInVCF = this.tmpDir + vcfName;
        String tmpVCF;
        Job preProcessVCF = getWorkflow().createBashJob("subset_VCF");
        Command cmd = preProcessVCF.getCommand();
        if (newInVCF.endsWith("gz")) {
            tmpVCF = newInVCF.replace(".vcf.gz", ".temp.vcf");
            cmd.addArgument("zcat " + inVCF + " >" + tmpVCF + ";\n");

        } else {
            tmpVCF = newInVCF.replace(".vcf", ".temp.vcf");
            cmd.addArgument("cp " + inVCF + " " + tmpVCF + ";\n");
        }
        cmd.addArgument(bedtools + "/bedtools intersect -header -a "
                + tmpVCF + " -b "
                + this.targetBedFile + " > "
                + tmpVCF.replace(".vcf", ".TGL.targ.vcf") + ";\n");
        String interTargVCF = tmpVCF.replace(".vcf", ".TGL.targ.vcf");
        cmd.addArgument(this.bgzip + " -c " + interTargVCF + " > "
                + interTargVCF + ".gz" + ";\n");
        cmd.addArgument(this.tabix + " -p vcf " + interTargVCF + ".gz");
        preProcessVCF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        preProcessVCF.setQueue(getOptionalProperty("queue", ""));
        hmap.put(interTargVCF + ".gz", preProcessVCF);
        return hmap;
    }

    private Job extractSampleNames(String inVCF) {
        Job extractSampleNames = getWorkflow().createBashJob("get_sample_ids");
        Command cmd = extractSampleNames.getCommand();
        cmd.addArgument(pathfixExport);
        cmd.addArgument("module load vcftools; \n");
        cmd.addArgument(vcftools + "/vcf-query -l " + inVCF + "> "
                + this.tmpDir + "sample_headers;\n");
        cmd.addArgument("cat " + this.tmpDir + "sample_headers"
                + " | grep -v \"GATK\" | tr \"\\n\" \",\" > "
                + this.tmpDir + "sample_names");
        extractSampleNames.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        extractSampleNames.setQueue(getOptionalProperty("queue", ""));
        return extractSampleNames;
    }
}
