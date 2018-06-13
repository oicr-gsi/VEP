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
    
    //vcf2maf
    private String vcf2mafpl;
    private String perlPath;
    private String vcf2mafPath;
    private String perl;
    
    //params
    private String hgBuild;
    private String species;
//    private String annotInfo = ",gnomAD,vcf,exact,0,AF_POPMAX,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS";
    private Integer hgvsShift = 1;
    private double vafFilter = 0.7;
    private String VEPpath;
    private String VEPdata;
//    private String retainInfo = "gnomAD_AF_POPMAX,gnomAD_AF_AFR,gnomAD_AF_AMR,gnomAD_AF_ASJ,gnomAD_AF_EAS,gnomAD_AF_FIN,gnomAD_AF_NFE,gnomAD_AF_OTH,gnomAD_AF_SAS";
    


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
            
            // vcf2maf
            vcf2mafpl = getProperty("VCF2MAF");
            perlPath = getProperty("PERL_PATH");
            vcf2mafPath=getProperty("VCF2MAF_PATH");
            perl=getProperty("TGL_PERL");

            // ref fasta
            refFasta = getProperty("ref_fasta");
            exacVCF = getProperty("ExAC_vcf");
            
            
            // VEP
            hgvsShift = Integer.parseInt(getOptionalProperty("hgvs_shift_flag", "1"));
            vafFilter = Double.parseDouble(getOptionalProperty("vaf_filter", "0.7"));
            species = getOptionalProperty("species", "homo_sapiens");
            hgBuild = getOptionalProperty("hg_version", "GRCh37");
            VEPpath = getProperty("VEP_PATH");
            VEPdata = getProperty("VEP_DATA");
            

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
//        this.outDir = this.outputFilenamePrefix + "_output/";
        String inVCF = getFiles().get("inVCF").getProvisionedPath();
        String tmpVCF = this.tmpDir + this.outputFilenamePrefix + ".temp.vcf";
        
        String mafFile = this.dataDir + this.outputFilenamePrefix + ".maf.txt";
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

        // Job pending to do: annotate OICR frequency
        
        Job vcf2MAF = runVcf2Maf(tmpVCF, mafFile);
        vcf2MAF.addParent(parentJob);
        parentJob = vcf2MAF;
            
        
        // Provision out maf.txt file
        SqwFile outMaf = createOutputFile(mafFile, TXT_METATYPE, this.manualOutput);
        outMaf.getAnnotations().put("MAF", "VEP");
        vcf2MAF.addFile(outMaf);
        
    }
          
    private Job runVcf2Maf(String inVCF, String outputMAF){
        Job runVCF2MAF = getWorkflow().createBashJob("vcf2maf");
        String PATHFIX = "# perl/5.22.2-tgl\n";
        PATHFIX =PATHFIX + "export LD_LIBRARY_PATH="+this.perlPath+"/lib:$LD_LIBRARY_PATH;\n";
        PATHFIX =PATHFIX + "export PERL5LIB="+this.perlPath+"/lib:$PERL5LIB\n";
        PATHFIX =PATHFIX + "export PATH="+this.perlPath+"/bin:$PATH\n";
        PATHFIX =PATHFIX + "# vep/92\n";
        PATHFIX =PATHFIX + "export PATH="+this.VEPpath+":$PATH\n";
        PATHFIX =PATHFIX + "export PATH="+this.VEPpath+"/htslib:$PATH\n";
        PATHFIX =PATHFIX + "export PATH="+this.VEPpath+"/samtools/bin:$PATH\n";
        PATHFIX =PATHFIX + "export PERL5LIB="+this.VEPpath+":$PATH\n";
        PATHFIX =PATHFIX + "export VEP_PATH="+ this.VEPpath + ";\n";
        PATHFIX =PATHFIX + "export VEP_DATA="+ this.VEPdata + ";\n";
        PATHFIX =PATHFIX + "# vcf2maf\n";
        PATHFIX =PATHFIX + "export PATH="+this.vcf2mafPath+":$PATH;\n";
        Command cmd = runVCF2MAF.getCommand();
//        cmd.addArgument("echo $MODULEPATH;");
//        cmd.addArgument("module use /.mounts/labs/PDE/Modules/modulefiles;");
//        cmd.addArgument("module load perl/5.22.2-tgl;");
//        cmd.addArgument("module load vep/92;");
//        cmd.addArgument("module load vcf2maf;");
//        cmd.addArgument("echo $MODULEPATH;");
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
        cmd.addArgument("--max-filter-ac 10");
        cmd.addArgument("--vep-path " + this.VEPpath);
        cmd.addArgument("--vep-data " + this.VEPdata);
//        cmd.addArgument("--retain-info " + this.retainInfo);
        cmd.addArgument("--min-hom-vaf "+ Double.toString(this.vafFilter));
        runVCF2MAF.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        runVCF2MAF.setQueue(getOptionalProperty("queue", ""));
        return runVCF2MAF;
    }
}