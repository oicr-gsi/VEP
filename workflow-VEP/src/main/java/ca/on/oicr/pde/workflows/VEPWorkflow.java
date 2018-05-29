package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
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
    private String outDir;

    // Input Data
    private String inputVCF;
    private String outputFilenamePrefix;
    private String normalSamplePrefix;
     

    //Tools
    private String samtools;
    private String java;
    private String gatk;
    private String bcftools;
    private String vcfscript;
    private String vcf2maf;
    private String vep;
    private String pythonpath;
    
//    private String bedtools;
    
    //VEP
    private String vepPath;
    private String vepData;
    
    // environment vars
    private String envVars;
    private String perl5lib;
    private String perlVersion = "5.10.1";
    
    
    //params
    private String hgBuild;
    private String species;
    private String annotInfo = ",gnomAD,vcf,exact,0,AF_POPMAX,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS";
    private Integer hgvsShift = 1;
    private double vafFilter = 0.7;
    private String retainInfo = "gnomAD_AF_POPMAX,gnomAD_AF_AFR,gnomAD_AF_AMR,gnomAD_AF_ASJ,gnomAD_AF_EAS,gnomAD_AF_FIN,gnomAD_AF_NFE,gnomAD_AF_OTH,gnomAD_AF_SAS";
    


    //Memory allocation
    private Integer VEPMem;



    //ref Data
    private String refFasta;
    private String exacVCF;
    private String gnomadVCF;
    
    private boolean manualOutput;
    private String queue;

    
//    // metatypes
    private String TXT_METATYPE=" application/txt-gz";
//    private String PDF_METATYPE="application/pdf";

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            inputVCF = getProperty("input_vcf_file");
            
            //
            outputFilenamePrefix = getProperty("external_identifier");
            normalSamplePrefix = getOptionalProperty("normal_sample_external_identifier", "NA");

            //tools
            vcf2maf = getProperty("vcf2maf");
            vep = getProperty("VEP");
            samtools = getProperty("samtools");

            // ref fasta
            refFasta = getProperty("ref_fasta");
            exacVCF = getProperty("ExAC_vcf");
            gnomadVCF = getProperty("GNOMAD_vcf");
            
            
            // VEP
            vepPath = getProperty("VEP_PATH");
            vepData = getProperty("VEP_DATA");
            hgvsShift = Integer.parseInt(getOptionalProperty("hgvs_shift_flag", "1"));
            vafFilter = Double.parseDouble(getOptionalProperty("vaf_filter", "0.7"));
            species = getOptionalProperty("species", "homo_sapiens");
            hgBuild = getOptionalProperty("hg_version", "GRCh37");
            
            // Environment vars
            perl5lib=getProperty("perl5lib");
            envVars = "export VEP_PATH="+this.vepPath+";";
            envVars = envVars + "export VEP_DATA="+vepData+";";
            envVars = envVars + "export PERL5LIB=$VEP_DATA:"+this.perl5lib+"; ";
            envVars = envVars + "export PATH=$VEP_PATH/htslib:$PATH; ";
            envVars = envVars + "export PATH="+this.samtools+":$PATH; ";

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
            file0.setType("application/vcf-gz");
            file0.setIsInput(true);  }
        else{
            file0.setSourcePath(inputVCF);
            file0.setType("application/vcf");
        }
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {   
        
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output/";
        String inVCF = getFiles().get("inVCF").getProvisionedPath();
        String tmpVCF = this.tmpDir + this.outputFilenamePrefix + ".temp.vcf";
        
        String annoGNOMADvcf = this.tmpDir + this.outputFilenamePrefix + "_gnomad" + ".vcf";
        String mafFile = this.outDir + this.outputFilenamePrefix + ".maf.txt";
        if (normalSamplePrefix == null || normalSamplePrefix == "NA"){
            this.normalSamplePrefix = this.outputFilenamePrefix;
        }
        // preprocess VCF 
        Job preprocessVCF = getWorkflow().createBashJob("preprocess_VCF");
        Command cmd = preprocessVCF.getCommand();
        if (inVCF.endsWith("gz")){
            cmd.addArgument("zcat " + inVCF + " >" + tmpVCF);
           
        }else {
            cmd.addArgument("cp "+ inVCF + " " + tmpVCF);
        }
         preprocessVCF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
         preprocessVCF.setQueue(getOptionalProperty("queue", ""));
         parentJob = preprocessVCF;
    
        Job annotateGNOMAD = annotateGnomad(tmpVCF, annoGNOMADvcf);
        annotateGNOMAD.addParent(parentJob);
        parentJob = annotateGNOMAD;
        
        Job vcf2MAF = runVcf2Maf(annoGNOMADvcf, mafFile);
        vcf2MAF.addParent(parentJob);
        parentJob = vcf2MAF;
            
        
        // Provision out maf.txt file
        SqwFile outMaf = createOutputFile(mafFile, TXT_METATYPE, this.manualOutput);
        outMaf.getAnnotations().put("MAF", "VEP83");
        vcf2MAF.addFile(outMaf);
        
    }
      
    private Job annotateGnomad(String inVCF, String gnomadVCF){
        Job annoGnomad = getWorkflow().createBashJob("annotate_gnomad");
        Command cmd = annoGnomad.getCommand();
        cmd.addArgument("module load perl/"+perlVersion + ";");
        cmd.addArgument(this.envVars);
        cmd.addArgument("perl "+this.vep);
        cmd.addArgument("--species "+ this.species);
        cmd.addArgument("--assembly "+ this.hgBuild);
        cmd.addArgument("--offline");
        cmd.addArgument("--no_progress");
        cmd.addArgument("--everything");
        cmd.addArgument("--shift_hgvs "+ Integer.toString(this.hgvsShift));
        cmd.addArgument("--check_exisiting");
        cmd.addArgument("--check_alleles");
        cmd.addArgument("--total_length");
        cmd.addArgument("--allele_number");
        cmd.addArgument("--no_escape");
        cmd.addArgument("--xref_refseq");
        cmd.addArgument("buffer_size "+"200");
        cmd.addArgument("--dir "+this.vepData);
        cmd.addArgument("--fasta "+this.refFasta);
        cmd.addArgument("--input_file "+ inVCF);
        cmd.addArgument("--force_overwrite");
        cmd.addArgument("--custom "+ this.gnomadVCF + this.annotInfo);
        cmd.addArgument("--vcf");
        cmd.addArgument("output_file "+ gnomadVCF);
        annoGnomad.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        annoGnomad.setQueue(getOptionalProperty("queue", ""));
        return annoGnomad;
    }
    
    private Job runVcf2Maf(String annoGNOMADvcf, String outputMAF){
        Job runVCF2MAF = getWorkflow().createBashJob("vcf2maf");
        Command cmd = runVCF2MAF.getCommand();
        cmd.addArgument("module load perl/"+perlVersion+ ";");
        cmd.addArgument(this.envVars);
        cmd.addArgument("perl "+this.vcf2maf);
        cmd.addArgument("--species "+ this.species);
        cmd.addArgument("--ncbi-build" + this.hgBuild);
        cmd.addArgument("--input-vcf " + annoGNOMADvcf);
        cmd.addArgument("--output-maf "+ outputMAF);
        cmd.addArgument("--tumor-id " + this.outputFilenamePrefix);
        cmd.addArgument("--normal-id " + this.normalSamplePrefix);
        cmd.addArgument("--vcf-tumor-id " + this.outputFilenamePrefix);
        cmd.addArgument("--vcf-normal-id " + this.normalSamplePrefix);
        cmd.addArgument("--vep-path "+ this.vepPath);
        cmd.addArgument("--vep-data "+ this.vepData);
        cmd.addArgument("--ref-fasta "+this.refFasta);
        cmd.addArgument("--filter-vcf "+this.exacVCF);
        cmd.addArgument("--max-filter-ac 10");
        cmd.addArgument("--retain-info " + this.retainInfo);
        cmd.addArgument("--min-hom-vaf "+ Double.toString(this.vafFilter));
        runVCF2MAF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        runVCF2MAF.setQueue(getOptionalProperty("queue", ""));
        return runVCF2MAF;
    }
}