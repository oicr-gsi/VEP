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
    


    //Memory allocation
    private Integer VEPMem;
    private String javaMem = "-Xmx16g";



    //ref Data
    private String refFasta;


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
            java = getProperty("java");
            vcfscript = getProperty("vcf_script").toString();

            // ref fasta
            refFasta = getProperty("ref_fasta");
            
            

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            VEPMem = Integer.parseInt(getProperty("picard_mem"));

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
        file0.setSourcePath(inputVCF);
        file0.setType("application/vcf-gz");
        file0.setIsInput(true);    
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output/";
        String inVCF = getFiles().get("inVCF").getProvisionedPath();
        
        String mafFile = this.outDir + this.outputFilenamePrefix + ".maf.txt";
//        String tumourSampleName = this.outputFilenamePrefix;
        if (normalSamplePrefix == null || normalSamplePrefix == "NA"){
            this.normalSamplePrefix = this.outputFilenamePrefix;
        }
        
        Job vepAnnotate = runVCFScript(inVCF, mafFile);
        
        // Provision out HS, HS2 and pdf
        SqwFile outMaf = createOutputFile(mafFile, TXT_METATYPE, this.manualOutput);
        outMaf.getAnnotations().put("MAF", "VEP83");
        vepAnnotate.addFile(outMaf);
        
    }
    
    
    private Job preProcVCF(String inVCF) {
        Job preProcess = getWorkflow().createBashJob("vcf_handling");
        Command cmd = preProcess.getCommand();
        if (inVCF.endsWith(".gz")){
        cmd.addArgument("zcat " + inVCF + " >" + this.outDir + this.outputFilenamePrefix + ".temp.vcf");
        }
        else{
            cmd.addArgument("ln -s "+inVCF + " " + this.outDir+this.outputFilenamePrefix + ".temp.vcf");
        }
        preProcess.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        preProcess.setQueue(getOptionalProperty("queue", ""));
        return preProcess;
    }  
    
    
    
    
    
    
    
    
    private Job runVCFScript(String inVCF, String mafFile) {
        Job vcfScript = getWorkflow().createBashJob("vcf_script_job");
        Command cmd = vcfScript.getCommand();
        cmd.addArgument(this.vcfscript);
        cmd.addArgument("INPUT="+ inVCF);
        cmd.addArgument("OUTPT="+mafFile);
        cmd.addArgument("EXTT="+this.outputFilenamePrefix);
        cmd.addArgument("EXTN="+this.normalSamplePrefix);
        vcfScript.setMaxMemory(Integer.toString(this.VEPMem * 1024));
        vcfScript.setQueue(getOptionalProperty("queue", ""));
        return vcfScript;
    }  
}