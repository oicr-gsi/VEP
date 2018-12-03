package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
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
    private String outputFilenamePrefix;
    private String normalSamplePrefix = null;
    
    //vcf2maf
    private String vcf2mafpl;
    private String perlPath;
    private String vcf2mafPath;
    private String perl;
    
    // preprocess
    String bgzip;
    String tabix;
    String tabixVersion;
    
    //params
    private String hgBuild;
    private String species;
    private double vafFilter = 0.7;
    private Integer bufferSize = 200;
    private Integer acFilter = 10;
    private String VEPpath;
    private String VEPdata;
    private String retainInfo;
    private String oncoKBPath;
    
    // additional params
    private String additionalParams;
    
    // path var
    String PATHFIX;


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
    private String TXT_GZ_METATYPE="application/txt-gz";
    private String TXT_METATYPE="plain/txt";
    private String VCF_METATYPE="application/vcf";
    private String VCF_GZ_METATYPE="application/vcf-gz";

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            inputVCF = getProperty("input_vcf_file");
            outputFilenamePrefix = getProperty("external_identifier");
            normalSamplePrefix = getOptionalProperty("matched_normal_name", "unmatched");
            
            
            // vcf2maf
            vcf2mafpl = getProperty("VCF2MAF");
            perlPath = getProperty("PERL_PATH");
            vcf2mafPath=getProperty("VCF2MAF_PATH");
            perl=getProperty("TGL_PERL");
            oncoKBPath = getProperty("ONCOKB");
            
            // preprocess
            tabixVersion = getProperty("tabix_version");
            bgzip = getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip";
            tabix = getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix";

            // ref fasta
            refFasta = getProperty("ref_fasta");
            // exac vcf
            exacVCF = getProperty("ExAC_vcf");
            // target bed file
            targetBedFile = getProperty("target_bed");
            //tglfreq file
            freqTextFile = getOptionalProperty("freq_file", null);
            
            
            // VEP
            vafFilter = Double.parseDouble(getOptionalProperty("vaf_filter", "0.7"));
            species = getOptionalProperty("species", "homo_sapiens");
            hgBuild = getOptionalProperty("hg_version", "GRCh37");
            VEPpath = getProperty("VEP_PATH");
            VEPdata = getProperty("VEP_DATA");
            bufferSize = Integer.parseInt(getOptionalProperty("buffer_size", "200"));
            acFilter = Integer.parseInt(getOptionalProperty("max_ac_filter", "10"));
            
            //additional params
            additionalParams = getOptionalProperty("additional_arg", null);
            
            //
            
            PATHFIX = "# perl/5.22.2-tgl\n";
            PATHFIX = PATHFIX + "export LD_LIBRARY_PATH="+this.perlPath+"/lib:$LD_LIBRARY_PATH;\n";
            PATHFIX = PATHFIX + "export PERL5LIB="+this.perlPath+"/lib:$PERL5LIB\n";
            PATHFIX = PATHFIX + "export PATH="+this.perlPath+"/bin:$PATH\n";
            PATHFIX = PATHFIX + "# vep/92\n";
            PATHFIX = PATHFIX + "export PATH="+this.VEPpath+":$PATH\n";
            PATHFIX = PATHFIX + "export PATH="+this.VEPpath+"/htslib:$PATH\n";
            PATHFIX = PATHFIX + "export PATH="+this.VEPpath+"/samtools/bin:$PATH\n";
            PATHFIX = PATHFIX + "export PERL5LIB="+this.VEPpath+":$PATH\n";
            PATHFIX = PATHFIX + "export VEP_PATH="+ this.VEPpath + ";\n";
            PATHFIX = PATHFIX + "export VEP_DATA="+ this.VEPdata + ";\n";
            PATHFIX = PATHFIX + "# vcf2maf\n";
            PATHFIX = PATHFIX + "export PATH="+this.vcf2mafPath+":$PATH;\n";
            PATHFIX = PATHFIX + "export PATH="+ this.oncoKBPath +":$PATH;\n";
            

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
        if (inputVCF.endsWith("gz")){
            file0.setSourcePath(inputVCF);
            file0.setType(VCF_GZ_METATYPE);
            file0.setIsInput(true);  }
        else{
            file0.setSourcePath(inputVCF);
            file0.setType(VCF_METATYPE);
        }
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {   
        
        Job parentJob = null;
        String inVCF = getFiles().get("inVCF").getProvisionedPath();
        String mafFile = this.dataDir + this.outputFilenamePrefix + ".maf.txt";
        
        // extract sample names first
        Job extractSampleNames = this.extractSampleNames(inVCF);
        parentJob = extractSampleNames;
        String sampleHeaderFile = this.tmpDir + "sample_header";
        try {
            String sampleName = this.sampleNames(sampleHeaderFile);
            if (sampleName.contains(",")){
                String[] sampleIDs = sampleName.split(",");
                for (String s : sampleIDs){
                    if (s.contains("_R")){
                        this.normalSamplePrefix = s;
                    } else {
                        this.outputFilenamePrefix = s;
                    }
                }
            } else {
                this.outputFilenamePrefix = sampleName;
                this.normalSamplePrefix = "unmatched";
                Job preprocessUnmatchedVCF = handleUnmatchedVCF(inVCF);
                preprocessUnmatchedVCF.addParent(parentJob);
                parentJob = preprocessUnmatchedVCF;
                inVCF = this.tmpDir + this.outputFilenamePrefix + ".unmatched.vcf.gz";
            }
        } catch (IOException ex) {
            Logger.getLogger(VEPWorkflow.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        // subset VCF 
        HashMap<String,Job> getTargetVCF = preProcessVCF(inVCF);
        String subsetVCF = getTargetVCF.keySet().toArray()[0].toString();
        Job preprocess = getTargetVCF.get(subsetVCF);
        preprocess.addParent(parentJob);
        parentJob = preprocess;
        
        // provision out subset VCF
        SqwFile targetVCF = createOutputFile(subsetVCF, VCF_METATYPE, this.manualOutput);
        targetVCF.getAnnotations().put("Target_VCF", "VEP");
        parentJob.addFile(targetVCF);
        
        // annotate frequency
        String intVCF;
        if (this.freqTextFile != null){
            Job tglFreq = TGLFreqAnnotation(subsetVCF);
            intVCF = subsetVCF.replace(".vcf",".tglfreq.vcf");
            this.retainInfo = "TGL_Freq";
            tglFreq.addParent(parentJob);
            parentJob = tglFreq;
        } else {
            intVCF = subsetVCF;
        }
        
        // run vcf to maf
        Job vcf2MAF = runVcf2Maf(intVCF, mafFile);
        vcf2MAF.addParent(parentJob);
        parentJob = vcf2MAF;
        
        // oncokb annotator
        Job oncoKBAnnotate = getWorkflow().createBashJob("oncokb_annotate");
        //cmd.addArgument(PATHFIX);
        oncoKBAnnotate.setCommand(PATHFIX + this.oncoKBPath + "/" + "Mafannotator.py -i " + mafFile + " -o " + mafFile.replace(".txt", ".oncoKB.txt"));
        oncoKBAnnotate.addParent(parentJob);
        parentJob = oncoKBAnnotate;
        parentJob.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        parentJob.setQueue(getOptionalProperty("queue", ""));
        
        // zip maf file
        String oncoKBMafFile = mafFile.replace(".txt", ".oncoKB.txt");
        Job zipMafFile = getWorkflow().createBashJob("zip_maf");
        zipMafFile.setCommand(bgzip + " -c " + oncoKBMafFile + " " + mafFile + ".gz");
        zipMafFile.addParent(parentJob);
        parentJob = zipMafFile;
        parentJob.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        parentJob.setQueue(getOptionalProperty("queue", ""));

        // Provision out maf.txt file
        SqwFile outMaf = createOutputFile(mafFile + ".gz", TXT_GZ_METATYPE, this.manualOutput);
        outMaf.getAnnotations().put("MAF", "VEP");
        parentJob.addFile(outMaf);
        
    }
          
    private Job runVcf2Maf(String inVCF, String outputMAF){
        Job runVCF2MAF = getWorkflow().createBashJob("vcf2maf");
        Command cmd = runVCF2MAF.getCommand();
        cmd.addArgument(PATHFIX);
        cmd.addArgument(this.perl + " " + this.vcf2mafpl);
        cmd.addArgument("--species "+ this.species);
        cmd.addArgument("--ncbi-build " + this.hgBuild);
        cmd.addArgument("--input-vcf " + inVCF);
        cmd.addArgument("--output-maf "+ outputMAF);
        cmd.addArgument("--tumor-id " + this.outputFilenamePrefix);
        cmd.addArgument("--normal-id " + this.normalSamplePrefix);
        cmd.addArgument("--vcf-tumor-id " + this.outputFilenamePrefix);
        cmd.addArgument("--vcf-normal-id " + this.normalSamplePrefix);
        cmd.addArgument("--ref-fasta "+this.refFasta);
        cmd.addArgument("--filter-vcf "+this.exacVCF);
        cmd.addArgument("--max-filter-ac " + this.acFilter);
        cmd.addArgument("--vep-path " + this.VEPpath);
        cmd.addArgument("--vep-data " + this.VEPdata);
        if (this.freqTextFile != null){
            cmd.addArgument("--retain-info " + this.retainInfo);
        }
        cmd.addArgument("--min-hom-vaf "+ Double.toString(this.vafFilter));
        cmd.addArgument("--buffer-size " + Integer.toString(this.bufferSize));
        if (this.additionalParams != null){
            cmd.addArgument(this.additionalParams);
        }
        runVCF2MAF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        runVCF2MAF.setQueue(getOptionalProperty("queue", ""));
        return runVCF2MAF;
    }
    
    // merge VCFs for mutect2 unmatched
    private Job handleUnmatchedVCF(String inVCF){
        String tempTumorVCF = this.tmpDir + this.outputFilenamePrefix + ".vcf";
        String tempMutect2VCF = this.tmpDir + this.outputFilenamePrefix + "unmatched.vcf";
        Job mergeMutect2VCF = getWorkflow().createBashJob("preprocess_unmatched");
        Command cmd = mergeMutect2VCF.getCommand();
        cmd.addArgument("sed -i \"s/QSS\\,Number\\=A/QSS\\,Number\\=\\./\" " + inVCF + ";\n");
        cmd.addArgument("echo -e " + this.outputFilenamePrefix + "\\n" + this.normalSamplePrefix + " > " + this.tmpDir + this.outputFilenamePrefix + "_header");
        cmd.addArgument("module load bcftools \n");
        cmd.addArgument("bcftools  merge " + inVCF + " " + inVCF + "--force-samples >" + tempTumorVCF + ";\n");
        cmd.addArgument("bcftools reheader -s " + this.tmpDir + this.outputFilenamePrefix + "_header " + tempTumorVCF + ">" + tempMutect2VCF + ";\n");
        cmd.addArgument(bgzip + " < " + tempMutect2VCF + " > " + tempMutect2VCF + ".gz" + ";\n");
        cmd.addArgument(tabix + " -p vcf " + tempMutect2VCF + ".gz");
        mergeMutect2VCF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        mergeMutect2VCF.setQueue(getOptionalProperty("queue", ""));
        return mergeMutect2VCF;
    }
    
    private Job TGLFreqAnnotation(String inVCF){
        Job annotateTGLFreq = getWorkflow().createBashJob("tgl_freq");
        Command cmd = annotateTGLFreq.getCommand();
        cmd.addArgument(PATHFIX);
        cmd.addArgument("module load bcftools; \n");
        cmd.addArgument("bcftools annotate -a " + this.freqTextFile);
        cmd.addArgument("-c CHROM,POS,REF,ALT,TGL_Freq");
        cmd.addArgument("-h <(echo '##INFO=<ID=TGL_Freq,Number=.,Type=Float,Description=\"Variant Frequency Among TGL Tumours (MuTect2 Artifact Detection)\">'");
        cmd.addArgument(inVCF + " | " + this.bgzip + " -c >" + inVCF.replace(".vcf", ".temp.vcf.gz") + ";\n");
        cmd.addArgument("echo \"Marking novel variants as TGL_Freq=0.0\"\n");
        cmd.addArgument("bcftools annotate -a " + this.freqTextFile);
        cmd.addArgument("-c CHROM,POS,REF,ALT,TGL_Freq");
        cmd.addArgument("-m \"-TGL_Freq=0.0\" ");
        cmd.addArgument(inVCF.replace(".vcf", ".temp.vcf.gz") + " | grep -v \"Sites not listed in OICR_Freq=0.0\" > " + inVCF.replace(".vcf",".tglfreq.vcf"));
        annotateTGLFreq.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        annotateTGLFreq.setQueue(getOptionalProperty("queue", ""));
        return annotateTGLFreq;
    }
    
    
    private HashMap<String,Job> preProcessVCF(String inVCF){
        HashMap<String,Job> hmap = new HashMap<String,Job>();
        int index = inVCF.lastIndexOf("/");
        String vcfName = inVCF.substring(index + 1);
        String newInVCF = this.tmpDir + vcfName;
        String tmpVCF;
        Job preProcessVCF = getWorkflow().createBashJob("subset_VCF");
        Command cmd = preProcessVCF.getCommand();
        cmd.addArgument("module load bedtools; \n");
        if (newInVCF.endsWith("gz")){
            tmpVCF = newInVCF.replace(".vcf.gz", ".temp.vcf");
            cmd.addArgument("zcat " + inVCF + " >" + tmpVCF + ";\n");
           
        }else {
            tmpVCF = newInVCF.replace(".vcf", ".temp.vcf");
            cmd.addArgument("cp "+ inVCF + " " + tmpVCF + ";\n");
        }
        cmd.addArgument("bedtools intersect -header -a " 
                + tmpVCF + " -b " 
                + this.targetBedFile + " > " 
                + tmpVCF.replace(".vcf", ".TGL.targ.vcf") + ";\n");
        preProcessVCF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        preProcessVCF.setQueue(getOptionalProperty("queue", ""));
        hmap.put(tmpVCF, preProcessVCF);
        return hmap;
    }
    
    private Job extractSampleNames(String inVCF){
        Job extractSampleNames = getWorkflow().createBashJob("get_sample_ids");
        Command cmd = extractSampleNames.getCommand();
        cmd.addArgument("module load vcftools;\n");
        cmd.addArgument("vcf-query -l " + inVCF  + "> " + this.tmpDir + "sample_headers;\n");
        cmd.addArgument("cat " + this.tmpDir + "sample_headers" + " | grep -v \"GATK\" | tr \"\\n\" \",\" > " + this.tmpDir + "sample_names");
        extractSampleNames.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        extractSampleNames.setQueue(getOptionalProperty("queue", ""));
        return extractSampleNames;
    }
    
    private String sampleNames(String sampleHeader) throws IOException{
        List<String> lines = Files.readAllLines(Paths.get(sampleHeader));
        String str = lines.toString();
        return str.substring(0, str.length() - 1);
    }
}