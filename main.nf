#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/cnvcall
========================================================================================
 nf-core/cnvcall Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/cnvcall
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/cnvcall --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
     --input                        Path to input TSV file on mapping

    Options:
      --genome                      Name of iGenomes reference
      --singleEnd                   Specifies that the input is single end reads

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
================================================================================
                               CHECKING REFERENCES AND PARAMETERS
================================================================================
*/
// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
params.bwaIndex = params.genome && params.fasta ? params.genomes[params.genome].bwaIndex ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fastaFai = params.genome && params.fasta ? params.genomes[params.genome].fastaFai ?: null : null

ch_fasta = params.fasta  ? Channel.value(file(params.fasta)) : "null"
ch_fastaFai = params.fastaFai  ? Channel.value(file(params.fastaFai)) : "null"
ch_bwaIndex = params.bwaIndex ? Channel.value(file(params.bwaIndex)) : bwaIndexes

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}


// Parse the TSV file and generate the input BAMs for all the Callers

tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") )) tsvPath = params.input

inputSample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    inputSample = extractBam(tsvFile)
}  else exit 1, 'Need to define tsvFile  see --help'

(bamManta, bamTIDDIT, bamGRIDSS, bamLumpy) = inputSample.into(4)

/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Fasta Ref']        = params.fasta
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch') {
  summary['AWS Region']     = params.awsregion
  summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-cnvcall-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/cnvcall Workflow Summary'
    section_href: 'https://github.com/nf-core/cnvcall'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

// Stage config files
ch_output_docs = Channel.fromPath("${baseDir}/docs/output.md")

// STEP MANTA

process Manta {
    label 'cpus_max'
    label 'memory_max'
    
    container('quay.io/biocontainers/manta:1.6.0--py27_0')
     
    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publishDirMode

    input:
    
        set idPatient,gender,status, idSample, file(bam), file(bai) from bamManta
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta


    script:
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    
    configManta.py \
        --bam ${bam} \
        --reference ${fasta} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    #Here we may pick only SV final as output 
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}

vcfManta = vcfManta.dump(tag:'Single Manta')

// STEP TIDDIT 

process TIDDIT {
    container('tiddit:v1')
    label 'cpus_max'
    label 'memory_max'
    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/TIDDIT", mode: params.publishDirMode

    publishDir params.outdir, mode: params.publishDirMode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        set idPatient,gender,status, idSample, file(bam), file(bai) from bamTIDDIT
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfTIDDIT
        set file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig") into tidditOut

    script:
    """
    TIDDIT.py --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}
    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf
    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf
    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz
    tabix TIDDIT_${idSample}.vcf.gz
    """
}

vcfTIDDIT = vcfTIDDIT.dump(tag:'TIDDIT')


process Lumpy{

    publishDir "${params.outdir}/VariantCalling/${idSample}/Lumpy", mode: params.publishDirMode 
    container('lumpysv:0.2.13' )

    label 'cpus_max'
    label 'memory_max'

    tag {idSample}

    input:
       set idPatient,gender,status, idSample, file(bam), file(bai) from bamLumpy

    output:
       set val("Lumpy"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfLumpy

    script:

    """
        samtools sort  -n -O SAM ${bam}| samblaster --ignoreUnmated --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        | samtools view -S -b - \
        > ${idSample}.bam

        lumpyexpress \
            -B ${idSample}.bam \
            -o ${idSample}_lumpysv.vcf

        # Sort by chr and position the vcf
            

        svtyper \
            -B ${bam} \
            -i ${idSample}_lumpysv.vcf    > ${idSample}_lumpysv.gt.vcf
        
        awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' ${idSample}_lumpysv.gt.vcf > ${idSample}_lumpysv.gt.sorted.vcf

        mv ${idSample}_lumpysv.gt.sorted.vcf  ${idSample}_lumpysv.gt.vcf

        bgzip ${idSample}_lumpysv.gt.vcf
        tabix ${idSample}_lumpysv.gt.vcf.gz

    """


}
vcfLumpy = vcfLumpy.dump(tag:'Lumpy')


// Define the blacklist from Assets

if (params.genome=='GRCh38') blacklist =  new File("$baseDir/assets/hg38-blacklist.v2.bed")
if (params.genome=='GRCh37') blacklist =  new File("$baseDir/assets/hg19-blacklist.v2.bed")

process gridss{

    publishDir "${params.outdir}/VariantCalling/${idSample}/GRIDSS", mode: params.publishDirMode 
    container('gridds:2.8.0' )

    label 'cpus_max'
    label 'memory_max'

    tag {idSample}

    input:
        set idPatient,gender,status, idSample, file(bam), file(bai) from bamGRIDSS
        file(bwaIndex) from ch_bwaIndex
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(blacklist) from blacklist

    output:
        set val("GRIDSS"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfGRIDSS

    script:
    """
        mkdir -p tmp
        
        gridss.sh \
            --reference ${fasta} \
            --jar \$GRIDSS_JAR \
            --output ${idSample}GRIDSS_CNV.vcf.gz \
            --assembly assembly.bam \
            --threads ${task.cpus} \
            --workingdir ./tmp/ \
            --maxcoverage 50000 \
            --label ${idSample} \
            --steps All \
            --blacklist ${blacklist} \
            --jvmheap ${task.memory.toGiga()}g \
            --picardoptions VALIDATION_STRINGENCY=LENIENT \
            ${bam}

    """

}


process combine_output {
    input:
       file(fasta) from ch_fasta
       set tiddit, idPatient_t, idSample_t,file(tiddit_vcf), file(tiddit_vcf_idx) from vcfTIDDIT
       set manta, idPatient_m, idSample_m, file(manta_vcf), file(manta_vcf_idx)   from vcfManta
       set lumpy, idPatient_l, idSample_l, file(lumpy_vcf), file(lumpy_vcf_idx)   from vcfLumpy
       set gridss, idPatient_l, idSample_g, file(gridss_vcf), file(gridss_vcf_idx) from vcfGRIDSS

    output:
       file('errorme') into bah_ch

    script:
    
    """
    echo 'error'
    exit 1
    """
     

}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/cnvcall] Successful: $workflow.runName"
    if (!workflow.success) {
      subject = "[nf-core/cnvcall] FAILED: $workflow.runName"
    }
    
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[nf-core/cnvcall]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/cnvcall]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/cnvcall v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

/*
================================================================================
                                 pipeline functions 
================================================================================
*/
// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [idSample, genderMap, statusMap, channel]
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, needed ${number} fields see --help for more information"
    return true
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
