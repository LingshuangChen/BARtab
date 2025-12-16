// Filters reads from BAM file that contain cell barcode (only contains R2 from CR)
// converts to fastq file
// https://kb.10xgenomics.com/hc/en-us/articles/360022448251-How-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
process BAM_TO_FASTQ {
    tag "${sample_id}"
    label "process_high_sc"
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam), path(index)

    output:
        tuple val(sample_id), path("${sample_id}_R2.fastq.gz"), emit: reads

    script:
        def umi_tag = params.pipeline == "splitpipe" ? "pN" : "UB"
        def unmapped = params.bam_unmapped ? "-f 4" : ""
        def contigs = params.bam_contigs ? params.bam_contigs.replace(",", " ") : ""
        // fix header if pipeline is cell ranger or similar (i.e. anything but splitpipe or saw)
        def pipeline = ["splitpipe", "saw"].contains(params.pipeline) ? params.pipeline : "starcode"
        def expr = "exists([CB]) && exists([${umi_tag}])"
        
        // Save the header lines
        // Filter alignments. Use LC_ALL=C to set C locale instead of UTF-8
        // can only filter for one tag in samtools
        // Combine header and body
        // convert BAM to fastq. CR output only contains R2
        // pipe everything to save time on IO


        // if both unmapped reads AND mapped reads have to be extracted, then they have to be concatenated.
        if(params.bam_contigs && params.bam_unmapped) {
        """
        cat <(samtools view -@ \$((${task.cpus} / 2)) -h $unmapped -e "$expr" $bam) \
            <(samtools view -@ \$((${task.cpus} / 2)) -e "$expr" $bam $contigs) |\
        if [[ "$pipeline" == "starcode" ]]; then
          samtools fastq -@ ${task.cpus} -T CB,UB |\
          sed -E "s/\tCB:Z:([A-Z]*)(-[1-9])?/_\\1/g;s/\tUB:Z:/_/g"
        elif [[ "$pipeline" == "splitpipe" ]]; then
          samtools fastq -@ ${task.cpus} |\
          sed 's/^\\([0-9]*_[0-9]*_[0-9]*\\)\\(.*\\)__\\([ACTG]*\\)__\\([^\t]*\\)/\\1\\2__\\3__\\4|\\1|\\3/g'
        else
          samtools fastq -@ ${task.cpus}
        fi | pigz -p ${task.cpus} > ${sample_id}_R2.fastq.gz
        """

        } else {

        """
        samtools view -h -@ ${task.cpus} $unmapped -e "$expr" $bam $contigs |\
        if [[ "$pipeline" == "starcode" ]]; then
          samtools fastq -@ ${task.cpus} -T CB,UB |\
          sed -E "s/\tCB:Z:([A-Z]*)(-[1-9])?/_\\1/g;s/\tUB:Z:/_/g"
        elif [[ "$pipeline" == "splitpipe" ]]; then
          samtools fastq -@ ${task.cpus} |\
          sed 's/^\\([0-9]*_[0-9]*_[0-9]*\\)\\(.*\\)__\\([ACTG]*\\)__\\([^\t]*\\)/\\1\\2__\\3__\\4|\\1|\\3/g'
        else
          samtools fastq -@ ${task.cpus}
        fi | pigz -p ${task.cpus} > ${sample_id}_R2.fastq.gz
        """
        }
}
