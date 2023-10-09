#!/bin/sh
source ~/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2018 Feb 9
######---------------------------------------------------------



######---------------------------------------------------------
######    Start the log
######---------------------------------------------------------

date | tee -a ~/scripthistory_gwas.log
echo 'Running imputation.sh with the following arguments:' \
		| tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------

######---------------------------------------------------------
######    Parse arguments and write them on the log
######---------------------------------------------------------

OPTIND=1
dataset=''
impprog=''
threads=8
memory=16000

while getopts "hd:i:t:m:" opt; do
	case "$opt" in
		h) echo "Performs individual and SNP quality checks on a dataset."
			echo "Usage: dataset_qc.sh -d [dataset] "
			echo "			-i [imputation software (impute/beagle/minimac)] "
			echo "			-t [threads (optional: default 8)] "
			echo "			-m [memory (optional: defaul 16000)]"
			sed -i '$d' ~/scripthistory_gwas.log
			sed -i '$d' ~/scripthistory_gwas.log
			exit 0 ;;
		d) dataset=`readlink -e $OPTARG.bed | sed 's/\.bed//'` 
			echo -e '\t -d '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		i) impprog=$OPTARG 
			echo -e '\t -i '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		t) threads=$OPTARG 
			echo -e '\t -t '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		m) memory=$OPTARG 
			echo -e '\t -m '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
	esac
done
echo '' | tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------

function pwait() {
	while [ $(jobs -p | wc -l) -ge $1 ]; do
		sleep 1
	done
}

######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_imp_${impprog}
cd ${dataset##*/}_imp_${impprog}

echo '\documentclass{article}

\usepackage[utf8]{inputenc}
\usepackage[a4paper, margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{authblk}
\usepackage{graphicx}
\usepackage{float}
\usepackage{listings}
\usepackage{tabularx}
\usepackage{array}
\usepackage{adjustbox}
\newcolumntype{L}{>{$}l<{$}}

\begin{document}
' > report.tex
echo '\documentclass{beamer}

\usepackage[utf8]{inputenc}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{float}
\usepackage{listings}
\usepackage{tabularx}
\usepackage{array}
\usepackage{adjustbox}
\newcolumntype{L}{>{$}l<{$}}

\usetheme{default}

\begin{document}
' > slides.tex


echo -e '\n'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞'\
		'🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞'\
		'🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜                                                             🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜                                                                   🙞🙔 🙞🙔'
echo -e '🙒🙜                                                                         🙞🙔'






######---------------------------------------------------------
######    Using SHAPEIT/IMPUTE2
######---------------------------------------------------------

if [ "$impprog" == "impute" ] ; then
	
	mkdir conform
	mkdir shapeit
	mkdir impute2
	
	cd conform
	
	for chr in $(seq 1 22) ; do
		
		echo "awk '{if(\$1==${chr}) print}' "${dataset}.bim | bash > ${dataset##*/}_chr${chr}.bim
		
		echo "zgrep 'rs' ${REFPATH}impute2/1000GP_Phase3/1000GP_Phase3_chr"${chr}".legend.gz \
			| awk '{print \$1}' | awk -F':' '{print ${chr},\$1,0,\$2,\$3,\$4}' " \
			| bash > ref_chr${chr}.bim
		
		python ${SCRIPTPATH}check-bims.py ref_chr${chr}.bim ${dataset##*/}_chr${chr}.bim \
			| tee bims_overlap_chr${chr}.log
		
		python ${SCRIPTPATH}bim-analyser.py ${dataset##*/}_chr${chr}.bim -ambiguous \
			> ${dataset##*/}_chr${chr}.bim.analysis.log
		
		cat ${dataset##*/}_chr${chr}.bim.ambiguous \
			ref_chr${chr}.${dataset##*/}_chr${chr}-problem.snps \
			<(grep DUP ${dataset##*/}_chr${chr}.bim) \
			<(grep -w 'I\|D' ${dataset##*/}_chr${chr}.bim) \
			> ${dataset##*/}_chr${chr}_removed.snps
		
		
		
		echo -e '\nStep ❶ for Chromosome '${chr}
		echo -e ${hru}
		${PLINK2PATH} --bfile ${dataset} \
						--allow-no-sex \
						--chr ${chr} \
						--flip ref_chr${chr}.${dataset##*/}_chr${chr}-needflip.snps \
						--exclude ${dataset##*/}_chr${chr}_removed.snps \
						--list-duplicate-vars \
						--out ${dataset##*/}_chr${chr}_impqc
		echo -e ${hrb}''
		echo -e '\n'
		
		
		
		cat <(awk '{print $4}{print $5}' ${dataset##*/}_chr${chr}_impqc.dupvar | grep -v 'rs') \
			<(grep .*rs.*rs ${dataset##*/}_chr${chr}_impqc.dupvar) \
			> ${dataset##*/}_chr${chr}_nonrsdups.snps
		
		cat ${dataset##*/}_chr${chr}.bim.ambiguous \
			ref_chr${chr}.${dataset##*/}_chr${chr}-problem.snps \
			${dataset##*/}_chr${chr}_nonrsdups.snps \
			<(grep DUP ${dataset##*/}_chr${chr}.bim) \
			> ${dataset##*/}_chr${chr}_removed.snps
		
		
		
		echo -e '\nStep ❷ for Chromosome '${chr}
		echo -e ${hru}
		${PLINK2PATH} --bfile ${dataset} \
						--allow-no-sex \
						--chr ${chr} \
						--flip ref_chr${chr}.${dataset##*/}_chr${chr}-needflip.snps \
						--exclude ${dataset##*/}_chr${chr}_removed.snps \
						--make-bed \
						--list-duplicate-vars \
						--out ${dataset##*/}_chr${chr}_forimp
		echo -e ${hrb}''
		echo -e '\n'
		
		
		
		echo -e '\nStep ❸ for Chromosome '${chr}
		echo -e ${hru}
		${SHAPEITPATH} --seed 12345 \
						-B ${dataset##*/}_chr${chr}_forimp \
						-O ${dataset##*/}_chr${chr}_shape \
						-L ${dataset##*/}_chr${chr}_shape.log \
						-T ${threads} \
						-M ${REFPATH}impute2/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
		echo -e ${hrb}''
		echo -e '\n'
		
		cat <(head -2 ${dataset##*/}_chr${chr}_shape.sample) \
			<(tail -n +3 ${dataset##*/}_chr${chr}_shape.sample \
				| awk '{print $1,$2,$3,$4,$5,$6,$7-1}' ) \
						> ${dataset##*/}_chr${chr}_shape.sample_fp
		
		gzip -f ${dataset##*/}_chr${chr}_shape.haps &
		
		
	done
	
	
	gzip conform/*.snps
	
	echo -e '\nStep ❹ '
	echo -e ${hru}
	for chr in $(seq 1 22) ; do
		
		for i in $(seq 1 5000000 $((($(awk '{print $4}' ${dataset##*/}_chr${chr}_forimp.bim \
			| tail -1)/5000000+1) * 5000000))) ; do
			 
			${IMPUTEPATH} -o_gz \
					-seed 1234567 \
					-use_prephased_g \
					-m ${REFPATH}impute2/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt \
					-h ${REFPATH}impute2/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz \
					-l ${REFPATH}impute2/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz \
					-known_haps_g ${dataset##*/}_chr${chr}_shape.haps.gz \
					-int $i $((i+5000000)) \
					-Ne 20000 \
					-k_hap 1000 \
					-buffer 1000 \
					-o ${dataset##*/}_chr${chr}_impute_$i &
			
			pwait ${threads}
		
		done
	
	done
	
	wait
	
	echo -e ${hrb}'\n'
	
	for i in ${dataset##*/}_chr*_impute_*.gz ; do
		
		if [ "$(zcat $i | wc -l)" -lt 1000 ] ; then 
			rm -rf $i
		fi
	done
	
	echo -e '\nStep ❺ '
	echo -e ${hru}
	for chr in $(seq 1 22) ; do
		
		for i in $(seq 1 5000000 $((($(awk '{print $4}' ${dataset##*/}_chr${chr}_forimp.bim \
			| tail -1)/5000000+1) * 5000000))) ; do
		
			if [ ! -f ${dataset##*/}_chr${chr}_impute_$i.gz ] ; then
				${IMPUTEPATH} -o_gz \
					-seed 1234567 \
					-use_prephased_g \
					-m ${REFPATH}impute2/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt \
					-h ${REFPATH}impute2/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz \
					-l ${REFPATH}impute2/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz \
					-known_haps_g ${dataset##*/}_chr${chr}_shape.haps.gz \
					-int $i $((i+5000000)) \
					-Ne 20000 \
					-buffer 1000 \
					-o ${dataset##*/}_chr${chr}_impute_$i &
			fi
			
			pwait ${threads}
		
		done
	
	done
	
	wait
	
	echo -e ${hrb}'\n'
	
	for i in grdiabetes_clean_chr*_impute_*1.gz ; do
		echo $i" "$(zcat $i | wc -l) >> imputed_snps_per_bin.list
	done
	
	for chr in $(seq 1 22) ; do
		ls impute2/${dataset##*/}_chr${chr}_impute*1.gz | sort -V | xargs zcat \
			| python ${SCRIPTPATH}impute_triple02na.py | gzip -c \
						> ${dataset##*/}_chr${chr}_impute_ALL.gz &
		pwait ${threads}
	done
	
	wait
	
	mkdir impute2_plink
	
	
	echo -e '\nStep ❻  - Convert from .gen to .bed'
	echo -e ${hru}
	
	for chr in $(seq 1 22) ; do
		
		(
		${PLINK2PATH} --allow-extra-chr \
						--allow-no-sex \
						--threads 1 \
						--memory ${memory} \
						--gen ${dataset##*/}_chr${chr}_impute_ALL.gz \
						--hard-call-threshold 0.1 \
						--make-bed \
						--out impute2_plink/${dataset##*/}_chr${chr}_imputed \
						--sample shapeit/${dataset##*/}_chr${chr}_shape.sample_fp
		
		sed -i "s/---/${chr}/g" impute2_plink/${dataset##*/}_chr${chr}_imputed.bim ) &
		
		pwait ${threads}
		
	done
	
	wait
	echo -e ${hrb}'\n'
	
	rm -rf ${dataset##*/}_chr*_impute_ALL.gz
	
	
	
	cd impute2_plink
	
	cat /dev/null > allchromosomes.mergelist
	
	for chr in $(seq 2 22) ; do
		echo ${dataset##*/}_chr${chr}_imputed.bed \
				${dataset##*/}_chr${chr}_imputed.bim \
				${dataset##*/}_chr${chr}_imputed.fam >> allchromosomes.mergelist
	done
	
	echo -e '\nStep ❼  - Merge all chromosomes together'
	echo -e ${hru}
	${PLINK2PATH} --bfile ${dataset##*/}_chr1_imputed \
					--merge-list allchromosomes.mergelist \
					--allow-no-sex \
					--threads ${threads} \
					--make-bed \
					--out ${dataset##*/}_chrALL_imputed
	echo -e ${hrb}'\n'
	
	for chr in $(seq 1 22) ; do
		rm -rf ${dataset##*/}_chr${chr}_imputed.bed \
				${dataset##*/}_chr${chr}_imputed.bim \
				${dataset##*/}_chr${chr}_imputed.fam
	done
	
	cat ../impute2/${dataset##*/}_chr*_impute_*_info | awk '{print $2,$7}' \
						| sort -V -u -k1,1 | grep -v 'rs_id' | gzip -fc > imputedsnps.info.gz &
	#cat ../${dataset##*/}_chr*_impute_*_info | awk '{print $2,1,2,$4,$5}' \
	#					> update.allele &
	cat ../impute2/${dataset##*/}_chr*_impute_*_info | awk '{if($7>0.1) print $2,$7}' \
						> over01.info &
	cat ../impute2/${dataset##*/}_chr*_impute_*_info | awk '{if($7>0.5) print $2,$7}' \
						> over05.info &
	cat ../impute2/${dataset##*/}_chr*_impute_*_info | awk '{if($7>0.8) print $2,$7}' \
						> over08.info &
	grep rs ${dataset##*/}_chrALL_imputed.bim | awk '{print $2}' > rsonly.snps
	wait
	
	echo -e '\nStep ❽ '
	echo -e ${hru}
	${PLINK2PATH} --bfile ${dataset##*/}_chrALL_imputed \
					--allow-no-sex \
					--freq \
					--missing \
					--out ${dataset##*/}_chrALL_imputed_check
	echo -e ${hrb}'\n'
	
	gzip -f ${dataset##*/}_chrALL_imputed_check.frq &
	gzip -f rsonly.snps &
	gzip -f ${dataset##*/}_chrALL_imputed_check.lmiss 
	wait
	
	zcat ${dataset##*/}_chrALL_imputed_check.frq | awk '{if($5<0.001) print $2,$5}' \
		| grep ':' > mafunder0001.snps
	zcat ${dataset##*/}_chrALL_imputed_check.frq | awk '{if($5<0.01) print $2,$5}' \
		| grep ':' > mafunder001.snps
	
	paste <( zcat imputedsnps.info.gz ) \
			<( zcat ${dataset##*/}_chrALL_imputed_check.frq.gz | tail -n +2 \
					| awk '{print $2,$5}' | sort -V -k1,1 ) \
			<( zcat ${dataset##*/}_chrALL_imputed_check.lmiss.gz | tail -n +2 \
					| awk '{print $2,$5}' | sort -V -k1,1 ) \
			| gzip -c > imputedstats.snps.gz
	python ${SCRIPTPATH}scatter_plot_simple.py imputedstats.snps.gz 1 3
	python ${SCRIPTPATH}scatter_plot_simple.py imputedstats.snps.gz 1 5
	
	echo -e '\nStep ❾ '
	echo -e ${hru}
	${PLINK2PATH} --bfile ${dataset##*/}_chrALL_imputed \
					--allow-no-sex \
					--extract over05.info \
					--exclude mafunder0001.snps \
					--threads ${threads} \
					--make-bed \
					--geno 0.05 \
					--out ${dataset##*/}_chrALL_imputed_info05_maf0001_geno005
	echo -e ${hrb}'\n'
	
	echo -e '\nStep ❾ '
	echo -e ${hru}
	${PLINK2PATH} --bfile ${dataset##*/}_chrALL_imputed \
					--allow-no-sex \
					--extract over05.info \
					--exclude mafunder001.snps \
					--threads ${threads} \
					--make-bed \
					--geno 0.02 \
					--out ${dataset##*/}_chrALL_imputed_info05_maf001_geno002
	echo -e ${hrb}'\n'
	
	gzip -f over01.info &
	gzip -f over05.info &
	gzip -f over08.info &
	gzip -f mafunder0001.snps &
	gzip -f mafunder001.snps
	wait
	
	cd ..
	
	#${GTOOLPATH} -G \
	#			--g ${dataset##*/}_chr${chr}_impute_ALL.gz \
	#			--s ${dataset##*/}_chr${chr}_shape.sample \
	#			--ped impute2_plink/${dataset##*/}_chr${chr}_impute_ALL.ped \
	#			--map impute2_plink/${dataset##*/}_chr${chr}_impute_ALL.map \
	#			--threshold 0.9 \
	#			--phenotype plink_pheno \
	#			--sex sex \
	#			--chr ${chr} & 
	#gzip -f gtool.log
	#Keep SNPs with MAF>0.005 and INFO>0.1
	#BEST GUESS GENOTYPES
	#genotype probs over 0.8, else NoCall
	#INFO Score (2*variance of dosages / ((1-frq)*frq)
	#split into noqc, subtle qc, strong qc
	#subtle qc call rate > .98
	#strong qc call rate > .99, MAF>.05
	
	echo '
	\newpage
	\section{Imputation QC of \texttt{\detokenize{'${dataset##*/}'}}}
		
	\begin{table}[H]
		\caption{}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{l|L|L|L|L|L}\toprule
			'`echo "CHR&BeforeCheck&AfterCheck&1000G(rs)&Imputed&I0.5M0.001G0.05&
				I0.5M0.01G0.02\\\\\\\\ \\\\hline" && \
				paste <(wc -l ${dataset##*/}_chr*[0-9].bim | sort -Vk2,2) \
				<(wc -l ${dataset##*/}*_forimp.bim | sort -Vk2,2) \
				<(wc -l conform/ref*.bim | sort -Vk2,2) \
				<(cat <(awk '{print $1}' impute2_plink/${dataset##*/}_chrALL_imputed.bim \
							| uniq -c ) \
						<( wc -l impute2_plink/${dataset##*/}_chrALL_imputed.bim ) ) \
				<(cat <(awk '{print $1}' \
					impute2_plink/${dataset##*/}_chrALL_imputed_info05_maf0001_geno005.bim \
					| uniq -c ) \
					<( wc -l \
					impute2_plink/${dataset##*/}_chrALL_imputed_info05_maf0001_geno005.bim ) ) \
				| awk '{print $6"&"$1"&"$3"&"$5"&"$7"&"$9"&"$11"\\\\\\\\ \\\\hline"}' \
				| sed -e 's/ref_//' -e 's/.bim//'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
			readlink -e impute2_plink/scatter-imputedstats.snps-13.png \
			| sed 's/\.png//' `'}.png}
		\caption{Scatter plot MAF(x) vs INFO(y)}
	\end{figure}
	
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
			readlink -e impute2_plink/scatter-imputedstats.snps-15.png \
			| sed 's/\.png//' `'}.png}
		\caption{Scatter plot GENO(x) vs INFO(y)}
	\end{figure}
	
	' >> report.tex
	
	echo '
	
	\section{Imputation QC of \texttt{\detokenize{'${dataset##*/}'}}}
	
	\begin{frame}{Imputation QC of \texttt{\detokenize{'${dataset##*/}'}}}
	\begin{table}[H]
		\caption{}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{l|L|L|L|L|L}\toprule
			'`echo "CHR&BeforeCheck&AfterCheck&1000G(rs)&Imputed&I0.5M0.001G0.05&
				I0.5M0.01G0.02\\\\\\\\ \\\\hline" && \
				paste <(wc -l ${dataset##*/}_chr*[0-9].bim | sort -Vk2,2) \
				<(wc -l ${dataset##*/}*_forimp.bim | sort -Vk2,2) \
				<(wc -l conform/ref*.bim | sort -Vk2,2) \
				<(cat <(awk '{print $1}' impute2_plink/${dataset##*/}_chrALL_imputed.bim \
							| uniq -c ) \
						<( wc -l impute2_plink/${dataset##*/}_chrALL_imputed.bim ) ) \
				<(cat <(awk '{print $1}' \
					impute2_plink/${dataset##*/}_chrALL_imputed_info05_maf0001_geno005.bim \
					| uniq -c ) \
					<( wc -l \
					\impute2_plink/${dataset##*/}_chrALL_imputed_info05_maf0001_geno005.bim ) ) \
				| awk '{print $6"&"$1"&"$3"&"$5"&"$7"&"$9"&"$11"\\\\\\\\ \\\\hline"}' \
				| sed -e 's/ref_//' -e 's/.bim//'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	\end{frame}
	
	\begin{frame}{Imputation QC of \texttt{\detokenize{'${dataset##*/}'}}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'\
			`readlink -e impute2_plink/scatter-imputedstats.snps-13.png \
			| sed 's/\.png//' `'}.png}
		\caption{Scatter plot MAF(x) vs INFO(y)}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Imputation QC of \texttt{\detokenize{'${dataset##*/}'}}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
			readlink -e impute2_plink/scatter-imputedstats.snps-15.png \
			| sed 's/\.png//' `'}.png}
		\caption{Scatter plot GENO(x) vs INFO(y)}
	\end{figure}
	\end{frame}
	
	
	' >> slides.tex
fi






######---------------------------------------------------------
######    Using BEAGLE
######---------------------------------------------------------



if [ "$impprog" == "beagle" ] ; then
	for i in $(seq 1 22) ; do
		${PLINK2PATH} --bfile ${dataset} \
						--allow-no-sex \
						--recode vcf-iid bgz \
						--chr ${chr} \
						--reference 
						--out ${dataset}_chr${chr}_forimp
	done
#		${BEAGLEPATH}conform-gt.24May16.cee.jar ${REFPATH}
fi






######---------------------------------------------------------
######    Using MiniMac
######---------------------------------------------------------



if [ "$impprog" == "minimac" ] ; then
	for chr in $(seq 1 22) ; do
		${PLINK2PATH} --bfile ${dataset} \
						--allow-no-sex \
						--chr ${chr} \
						--make-bed \
						--out ${dataset}_chr${chr}_forimp
	done
		${SHAPEITPATH} --seed 12345 \
						-B ${dataset}_chr${chr}_forimp \
						--output-max outname.haps \
						--output-log outname.shape.log \
						--input-map genetic_map_file
		${SHAPEITPATH} -convert \
						--input-haps \
						--output-vcf .vcf
fi






######---------------------------------------------------------
######     End of analysis and LaTeX report
######---------------------------------------------------------



echo -e '\n\n'
echo -e ${hru}
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'

echo '\end{document}' >> report.tex
pdflatex -interaction=batchmode report.tex
rm -rf report.aux report.log

echo '\end{document}' >> slides.tex
pdflatex -interaction=batchmode slides.tex
rm -rf slides.aux slides.log slides.nav slides.out slides.snm slides.toc
echo -e ${hrb}'\n'

echo -e '🙐🙝                                                                         🙟🙖'
echo -e '🙐🙝 🙐🙝                                                                   🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝                                                             🙟🙖 🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟'\
		'🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟'\
		'🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖'
echo -e '\n'
cd ..
