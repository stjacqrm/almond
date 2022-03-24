#!/usr/bin/env nextflow

//Description: Analyzing SC2 Competencies
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov


//setup channel to read in the fasta files
Channel
    .fromPath("${params.reference}")
    .ifEmpty { exit 1, "Cannot find the reference fasta file in ${params.reference}" }
    .set { reference }

if (params.blindComp=="BlindC"){
    Channel
        .fromPath("$baseDir/bam_files/BlindC-102221SK_consensus.fasta")
        .set { blindComp }

    Channel
        .fromPath("$baseDir/bam_files/BlindC-102221SK.sorted.bam")
        .set { blind_bam }
} else {
  if (params.blindComp=="BlindD"){
      Channel
          .fromPath("$baseDir/bam_files/BlindD-102221SK_consensus.fasta")
          .set { blindComp }

      Channel
          .fromPath("$baseDir/bam_files/BlindD-102221SK.sorted.bam")
          .set { blind_bam }
}else{
  if (params.blindComp=="BlindE"){
      Channel
          .fromPath("$baseDir/bam_files/BlindE-102221SK_consensus.fasta")
          .set { blindComp }

      Channel
          .fromPath("$baseDir/bam_files/BlindE-102221SK.sorted.bam")
          .set { blind_bam }
}}}


Channel
    .fromPath("${params.ptSample}")
    .ifEmpty { exit 1, "Cannot find the blindComp fasta file in ${params.ptSample}" }
    .set { ptSample }

Channel
    .fromPath("${params.bam_files}")
    .ifEmpty { exit 1, "Cannot find the sample bam files in ${params.bam_files}" }
    .set { bam_files }

Channel
    .fromPath("$baseDir/report/report_template.Rmd")
    .set { report }

Channel
    .fromPath("$baseDir/bash/create_almond_report.sh")
    .set { bash }

Channel
    .fromPath("$baseDir/bam_files/artic_V3_nCoV-2019.bed")
    .set { bed_files }


//Step 1: collect the SC2 assembled_genomes
process collect_fasta {
  publishDir "${params.outdir}/logs/multi_fasta/", mode: 'copy'


  input:
  file( reference ) from reference
  file( blindComp ) from blindComp
  file( ptSample ) from ptSample

  output:
  file("multi.fasta") into mafft, pangolin

  script:
  """
  cat ${reference} ${blindComp} ${ptSample} > multi.fasta;
  """

}

//Step 2: generate a multi alignment with mafft
process alignment {
  publishDir "${params.outdir}/results/mafft", mode: 'copy', pattern:'*'

  input:
  file(assembly) from mafft

  output:
  file("multi.cleaned.alignment") into multi_alignment
  //file("output.vcf") into output_vcf

  script:
  """
  mafft ${assembly} > multi.alignment; linenum=\$(grep -n "NC_045512.2" multi.alignment | cut -d':' -f1 | tail -1); sed -n "\$((linenum))"',\$p' multi.alignment > multi.cleaned.alignment
  """

}

process snp_sites {
  publishDir "${params.outdir}/results/snp_sites", mode: 'copy', pattern:'*'

  input:
  file(multi_alignment) from multi_alignment

  output:
  //file("multi.cleaned.alignment") into multi_alignment
  file("output.csv") into output_vcf,render_vcf

  script:
  """
  snp-sites -v -o output.vcf ${multi_alignment}
  sed -i -e '1,3d' output.vcf
  sed 's/\t/,/g' output.vcf > output.csv
  """

}

//Step 5: pangolin for lineages
process pangolin {
  publishDir "${params.outdir}/results/pangolin_results", mode: 'copy' , pattern:"*_report.csv"

  input:
  file(assembly) from pangolin

  output:
  file("*report.csv") into lineage_results,render_lineage
  file("*complete_lineage_version.csv") into pangolin_version

  script:
  """
  pangolin ${assembly} --outfile lineage_report.csv
  pangolin --version > pango_version.csv
  echo "version" >> complete_lineage_version.csv
  cat pango_version.csv >> complete_lineage_version.csv
  """
}

//Step 5: ampliconstats
process ampliconstats {
  publishDir "${params.outdir}/results/ampliconstats", mode: 'copy' , pattern:"*"

  input:
  file(pt_bam) from bam_files
  file(blind_bam) from blind_bam
  file(bed) from bed_files

  output:
  file("pt-sample-combined-reads.png") into pt_amp_stats
  file("blind-sample-combined-reads.png") into blind_amp_stats

  script:
  """
  samtools ampliconstats ${bed} ${pt_bam} 2>> errfile > pt_ampliconstats.txt
  plot-ampliconstats -size 1200,900 pt-sample pt_ampliconstats.txt

  samtools ampliconstats ${bed} ${blind_bam} 2>> errfile > blind_ampliconstats.txt
  plot-ampliconstats -size 1200,900 blind-sample blind_ampliconstats.txt
  """
}

//Step 7: generate the pdf report
process render{
  publishDir "${params.outdir}/", mode: 'copy', pattern:'Almond-report.pdf'
  //publishDir "${params.outdir}/", mode: 'copy', pattern:"berrywood*"
  //publishDir "${params.outdir}/results/berrywood_results/csvs", mode: 'copy' , pattern:"*csv"
  beforeScript 'ulimit -s unlimited'

  input:
  file(lineage) from render_lineage
  file(csv) from render_vcf
  file(pt_png) from pt_amp_stats
  file(blind_png) from blind_amp_stats
  file(rmd) from report
  file(bash) from bash

  output:
  //file('*.csv')
  file "Almond-report.pdf"

  shell:
  """
  cp ${rmd} ./report_template.Rnw
  cp ${bash} ./create_report.sh
  chmod +x create_report.sh
  bash create_report.sh -p "Almond-report" -t "${params.title}" -T report_template.Rnw -o . -a ${lineage} -s ${csv} -P ${pt_png} -b ${blind_png}
  """

}
