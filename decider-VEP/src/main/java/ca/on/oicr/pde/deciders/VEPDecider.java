package ca.on.oicr.pde.deciders;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author prath@oicr.on.ca
 */
public class VEPDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private String[] allowedTemplateTypes = {"EX", "WT"};
    private String templateType;
    private String queue = "";
//    private String outputFileNamePrefix;

    private String[] allowedExtensions = new String[]{".snv.indel.vcf.gz", ".tumor_only.vcf.gz"};
    private String extensions;

    // VEP
    private String exacVCF = "/oicr/local/analysis/sw/vep/vep92/.cache/Plugins/ExAC.r0.3.sites.minus_somatic.vcf.gz";
    private String hgBuild = "GRCh37";
    private String species = "homo_sapiens";
    private String vafFilter = "0.7"; //for homozygous calls
    private String bufferSize = "200";
    private String maxACFilter = "10";
    private String additionalArgs = null;

    // onco KB
    private String oncoKBpath = "/.mounts/labs/PDE/Modules/sw/oncokb-annotator";

    // memory
    private String vepMem = "30";

    private String targetBed = "/.mounts/labs/PDE/data/reference/targets/ensembl_v6_ccds_exons_onebed_intersectRegions_pad50.bed";
    private final String refGenome = "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19/Genomic/references/fasta/hg19.fa";
    private String freqDB = null;

    // metatype info
    private final static String VCF_METATYPE = "application/vcf";
    private final static String VCF_GZ_METATYPE1 = "application/vcf-4-gzip";
    private final static String VCF_GZ_METATYPE2 = "application/vcf-gz";
    private final static ArrayList<String> METATYPES = new ArrayList<String>(Arrays.asList(VCF_METATYPE, VCF_GZ_METATYPE1, VCF_GZ_METATYPE2));

    public VEPDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type. Default: " + String.join(",", this.allowedTemplateTypes)).withOptionalArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("target-bed", "Optional parameter for VEP workflow: Specify the path to interval bed file." + "Default is " + this.targetBed).withOptionalArg();
        parser.accepts("ref-fasta", "Optional parameter for VEP workflow: Specify the path to reference human genome fasta. Default: " + this.refGenome).withOptionalArg();
        parser.accepts("tgl-freq-file", "Optional parameter for VEP workflow: Specify the path to the file containing frequency information. Default: null").withOptionalArg();
        parser.accepts("exac-vcf", "Optional parameter for VEP workflow: Specify the path to the EXAC VCF. Default: " + this.exacVCF).withOptionalArg();
        parser.accepts("hg-build", "Optional parameter for VEP workflow: Specify the build of human genome reference. Default: " + this.hgBuild).withOptionalArg();
        parser.accepts("species", "Optional parameter for VEP workflow: Specify the species name. Default: " + this.species).withOptionalArg();
        parser.accepts("hom-vaf-filter", "Optional parameter for VEP workflow: Specify the minimum vaf for homozygous calls. Default: " + this.vafFilter).withOptionalArg();
        parser.accepts("buffer-size", "Optional parameter for VEP workflow: Specify the buffer size. Default: " + this.bufferSize).withOptionalArg();
        parser.accepts("max-ac-filter", "Optional parameter for VEP workflow: Specify the max AC filter. Default: " + this.maxACFilter).withOptionalArg();
        parser.accepts("additional-args", "Optional: Include additional arguments and paramenters to run VEP. Default: " + this.additionalArgs).withOptionalArg();
        parser.accepts("allow-extensions", "Optional:comma-separated VCF file extensions to annotate. Default " + String.join(",", this.allowedExtensions)).withOptionalArg();
        parser.accepts("vep-mem", "Optional: specify memory for VEP. Default " + this.vepMem).withOptionalArg();
        parser.accepts("oncokb-path", "Optional: Specify the path to the folder containing the script oncoKBAnnotator.py. Default. " + this.oncoKBpath).withOptionalArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(METATYPES);
        this.setHeadersToGroupBy(Arrays.asList(Header.FILE_SWA));

        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            this.templateType = options.valueOf("template-type").toString();
            this.allowedTemplateTypes = this.templateType.split(",");
        } else {
            Log.warn("Please re-run with --template-type <template type>. Supported template type(s) are EX,WT");
        }

        if (this.options.has("tgl-freq-file")) {
            this.freqDB = options.valueOf("tgl-freq-file").toString();
        }

        if (this.options.has("exac-vcf")) {
            this.exacVCF = options.valueOf("exac-vcf").toString();
        }

        if (this.options.has("hg-build")) {
            this.hgBuild = options.valueOf("hg-build").toString();
        }

        if (this.options.has("species")) {
            this.species = options.valueOf("species").toString();
        }

        if (this.options.has("hom-vaf-filter")) {
            this.vafFilter = options.valueOf("hom-vaf-filter").toString();
        }

        if (this.options.has("buffer-size")) {
            this.bufferSize = options.valueOf("buffer-size").toString();
        }

        if (this.options.has("max-ac-filter")) {
            this.maxACFilter = options.valueOf("max-ac-filter").toString();
        }

        if (this.options.has("additional-args")) {
            this.additionalArgs = options.valueOf("additional-args").toString();
        }

        if (this.options.has("target-bed")) {
            this.targetBed = options.valueOf("target-bed").toString();
        }

        if (this.options.has("allow-extensions")) {
            this.extensions = options.valueOf("allow-extensions").toString();
            this.allowedExtensions = this.extensions.split(",");
        }
        if (this.options.has("vep-mem")) {
            this.vepMem = options.valueOf("vep-mem").toString();
        }
        if (this.options.has("oncokb-path")) {
            this.oncoKBpath = options.valueOf("oncokb-path").toString();
        }
        return rv;
    }

    /**
     * Final check
     *
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     *
     * @return
     */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        Log.info(filePaths);
        if (filePaths.length != 1) {
            Log.info(filePaths.length);
            Log.info(filePaths);
            Log.error("Workflow runs can be scheduled for exactly 1 VCF");
        }

        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
        String resequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String target_bed = null;

        if (null == currentTissueType) {
            return false; // we need only those which have their tissue type set
        }

        if (currentTissueType.equals("R")) {
            return false;
        }

        if (!Arrays.asList(this.allowedTemplateTypes).contains(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Log.debug("Number of files from file provenance = " + vals.size());

        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            boolean fileExtensionOK = false;
            boolean templateTypeCheck = false;
            String currentTtype = currentRV.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            if (Arrays.asList(this.allowedTemplateTypes).contains(currentTtype)) {
                templateTypeCheck = true;
            }

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                String filePath = currentRV.getFiles().get(f).getFilePath();
                try {
                    if (METATYPES.contains(currentRV.getFiles().get(f).getMetaType())) {
                        metatypeOK = true;
                    }
                    if (this.identifyFilePath(filePath)) {
                        fileExtensionOK = true;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }

            if (!templateTypeCheck) {
                continue;
            }

            if (!metatypeOK) {
                continue; // Go to the next value
            }

            if (!fileExtensionOK) {
                continue;
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);

            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String vcfPath = commaSeparatedFilePaths.split(",")[0];
        String outputFileNamePrefix = getOutputFileNamePrefix(vcfPath);

        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        iniFileMap.put("input_vcf_file", vcfPath);
        iniFileMap.put("output_filename_prefix", outputFileNamePrefix);
        iniFileMap.put("target_bed", this.targetBed);
        if (!this.queue.isEmpty()) {
            iniFileMap.put("queue", this.queue);
        }
        iniFileMap.put("ref_fasta", this.refGenome);
        iniFileMap.put("exac_vcf", this.exacVCF);
        iniFileMap.put("vaf_filter", this.vafFilter);
        iniFileMap.put("buffer_size", this.bufferSize);
        iniFileMap.put("max_ac_filter", this.maxACFilter);
        iniFileMap.put("species", this.species);
        iniFileMap.put("hg_version", this.hgBuild);
        if (this.additionalArgs != null) {
            iniFileMap.put("additional_args", this.additionalArgs);
        }
        if (this.freqDB != null) {
            iniFileMap.put("freq_file", this.freqDB);
        }
        iniFileMap.put("vep_mem", this.vepMem);
        iniFileMap.put("oncokb", this.oncoKBpath);
        return iniFileMap;
    }

    private boolean identifyFilePath(String filePath) {
        return Arrays.stream(allowedExtensions).anyMatch(entry -> filePath.endsWith(entry));
    }

    private String getOutputFileNamePrefix(String inVCF) {
        String baseName = FilenameUtils.getBaseName(inVCF);
        String[] bNames = baseName.split("\\.");
        String sampleName = bNames[0];
        return sampleName;
    }

    public static void main(String args[]) {
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(VEPDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));
    }

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String tubeID = null;

        private String extName = null;
        private String groupID = null;
        private String groupDescription = null;
        private String rootSampleName = null;
        private String workflowName = null;

        public String getRootSampleName() {
            return rootSampleName;
        }

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            workflowName = rv.getAttribute(Header.WORKFLOW_NAME.getTitle());
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + workflowName;
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            extName = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_external_name");
            rootSampleName = rv.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
            if (null == extName || extName.isEmpty()) {
                extName = rootSampleName;
            }
            groupID = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == groupDescription || groupDescription.isEmpty()) {
                groupDescription = "NA";
            }
            StringBuilder gba = new StringBuilder(fa.getDonor());
            gba.append(":").append(fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));

            String trs = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
            if (null != trs && !trs.isEmpty()) {
                gba.append(":").append(trs);
            }

            groupByAttribute = gba.toString() + ":" + extName + ":" + groupID + ":" + workflowName;
            path = rv.getFiles().get(0).getFilePath() + "";

        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getTissueType() {
            return tissueType;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public String getTubeID() {
            return tubeID;
        }

        public String getExtName() {
            return extName;
        }

        public String getGroupID() {
            return groupID;
        }

        public String getGroupDescription() {
            return groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }

    public static List<String> detectDuplicates(String commaSeparatedFilePaths) {

        String[] filePaths = commaSeparatedFilePaths.split(",");
        List<String> list = new ArrayList<String>();
        List<String> checker = new ArrayList<String>();

        for (String path : filePaths) {
            String baseName = makeBasename(path, ".bam");

            if (checker.contains(baseName) && !list.contains(path)) {
                list.add(path);
            } else {
                checker.add(baseName);
            }
        }

        return list.isEmpty() ? null : list;

    }

    /**
     * Utility function
     *
     * @param path
     * @param extension
     *
     * @return
     */
    public static String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }
}
