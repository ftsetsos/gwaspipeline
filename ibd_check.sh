#!/bin/sh
source ~/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2018 Jan 2
######---------------------------------------------------------



######---------------------------------------------------------
######    Start the log
######---------------------------------------------------------

date | tee -a ~/scripthistory_gwas.log
echo 'Running ibd_check.sh with the following arguments:' \
		| tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------

######---------------------------------------------------------
######    Parse arguments and write them on the log
######---------------------------------------------------------

OPTIND=1
dataset=''
outliers=${BLANKFILE}
threads=8
memory=16000

while getopts "hd:o:t:m:" opt; do
	case "$opt" in
		h) echo "Performs IBD check to identify relatedness. Prunes it first."
			echo "Usage: ibd_check.sh -d [dataset] "
			echo "			-o [outliers file (optional)] "
			echo "			-t [threads (optional: default 8)] "
			echo "			-m [memory (optional: default 16000)]"
			sed -i '$d' ~/scripthistory_gwas.log
			sed -i '$d' ~/scripthistory_gwas.log
			exit 0 ;;
		d) dataset=`readlink -e $OPTARG.bed | sed 's/\.bed//'` 
			echo -e '\t -d '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		o) outliers=`readlink -e $OPTARG` 
			echo -e '\t -o '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		t) threads=$OPTARG 
			echo -e '\t -t '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		m) memory=$OPTARG 
			echo -e '\t -m '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
	esac
done
echo '' | tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------


######---------------------------------------------------------
######    Parameters - TAKE SPECIAL CARE IF YOU CHANGE THESE
maf=0.05
hwe=0.001
snpmissing=0.02

mhcchrom=6
mhcstart=25000000
mhcend=35000000

chr8invchrom=8
chr8invstart=7000000
chr8invstop=13000000

chr17invchrom=17
chr17invstart=40900000
chr17invstop=44900000
######---------------------------------------------------------



######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_ibd_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_ibd_$(echo ${outliers##*/} | cut -d'.' -f1)

date >> runtime.info

echo -e $latexdoc > report.tex
echo -e $latexbeamer > slides.tex

commandvar="awk '{if((\$1==$mhcchrom && \$4>$mhcstart && \$4<$mhcend) || \
					(\$1==$chr8invchrom && \$4>$chr8invstart && \$4<$chr8invstop) || \
					(\$1==$chr17invchrom && \$4>$chr17invstart && \$4<$chr17invstop)) \
						print \$2, \$1, \$4}' ${dataset}.bim > mhc817.snps"
eval $commandvar
awk '{if($1>22 || $1==0 ) print}' ${dataset}.bim > nonautosomal.snps
cat mhc817.snps nonautosomal.snps > mhc817nonauto.snps

echo -e $startheader

echo -e '\n\nStep ❶'
echo -e ${hru}
	${PLINK2PATH} --bfile ${dataset} \
					--geno ${snpmissing} \
					--maf ${maf} \
					--hwe ${hwe} \
					--exclude mhc817nonauto.snps \
					--indep-pairwise 200 100 0.2 \
					--remove ${outliers} \
					--out ibd1 \
					--memory ${memory}
echo -e ${hrb}'\n'

############   If there are still more than 150K SNPs in the dataset,
############   then prune again using the same settings.

if [ "$(wc -l ibd1.prune.in | awk '{print $1}')" -gt 150000 ] ; then 
	
	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --bfile ${dataset} \
						--extract ibd1.prune.in \
						--indep-pairwise 200 100 0.2 \
						--remove ${outliers} \
						--out ibd2 \
						--memory ${memory}
	echo -e ${hrb}'\n'
	
	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --bfile ${dataset} \
						--extract ibd2.prune.in \
						--make-bed \
						--remove ${outliers} \
						--out ibd \
						--memory ${memory}
	echo -e ${hrb}'\n'
	
else
	
	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --bfile ${dataset} \
						--extract ibd1.prune.in \
						--make-bed \
						--remove ${outliers} \
						--out ibd \
						--memory ${memory}
	echo -e ${hrb}'\n'
	
fi

echo -e '\nStep ❸'
echo -e ${hru}
	${PLINK2PATH} --bfile ibd \
					--nonfounders \
					--freq \
					--remove ${outliers} \
					--out ibdfreq \
					--memory ${memory}
echo -e ${hrb}'\n'

echo -e '\nStep ❹'
echo -e ${hru}
	${PLINK2PATH} --bfile ibd \
					--genome full unbounded \
					--min 0.05 \
					--remove ${outliers} \
					--threads ${threads} \
					--out ibdgenome \
					--memory ${memory}
echo -e ${hrb}'\n'

############   Just for check, keep a list of couples related >0.1
############...Normal threshold is 0.1875

awk '{if($10>0.1) print}' ibdgenome.genome > ibdcheck01.rel
awk '{if($10>0.1875) print}' ibdcheck01.rel > ibdcheck02.rel
awk '{print $3,$4}' ibdcheck02.rel > ibdout.ind

awk '{if($10>0.95) print}' ibdcheck01.rel > ibdcheck_dup.rel        #### See the duplicates
awk '{if(($10>0.4)&&($10<0.6)) print}' ibdcheck01.rel \
										> ibdcheck_50.rel    ####See the siblings, parents
awk '{if($8>0.9) print}' ibdcheck_50.rel > ibdcheck_parents.rel
awk '{if($8<0.9) print}' ibdcheck_50.rel > ibdcheck_siblings.rel

cat <(awk '{print $1,$2}' ibdcheck01.rel) \
	<(awk '{print $3,$4}' ibdcheck01.rel) | sort | uniq -c \
			| sort -gk1 | awk '{print $2,$3,$1}' > most_pairings01_double.ind

awk '{if($3>1) print}' most_pairings01_double.ind > top01.ind

cat <(awk '{print $1,$2}' ibdcheck02.rel) \
	<(awk '{print $3,$4}' ibdcheck02.rel) | sort | uniq -c \
			| sort -gk1 | awk '{print $2,$3,$1}' > most_pairings02_double.ind

awk '{print $1,$2}' ibdcheck02.rel | sort | uniq -c | sort -gk1 \
	| awk '{print $2,$3,$1}' > most_pairings02.ind

awk '{print $1,$2}' ibdcheck01.rel | sort | uniq -c | sort -gk1 \
	| awk '{print $2,$3,$1}' > most_pairings01.ind

echo `grep -v 'HOMHOM' ibdcheck02.rel | wc -l `' pairs with over 18.75% relatedness' \
												> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_dup.rel |  wc -l `' duplicate pairs' >> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_parents.rel |  wc -l `' parent pairs' >> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_siblings.rel |  wc -l `' sibling pairs' >> ibd_check_summary.txt

echo -e '\n\n'
echo -e ${hru}
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}histogram.py ibdgenome.genome -col=9 -header -binsize=0.01 -log
python ${SCRIPTPATH}ibd_scatter_plot.py ibdcheck01.rel 6:7
python ${SCRIPTPATH}ibd_scatter_plot.py ibdcheck01.rel 7:8
echo -e ${hrb}''


############   Create a new fam file using the newly acquired information
############   

awk '{print $1,"Relatedness/"$1"-"$3"/("$7","$8","$9","$10")\n"$3,"Relatedness/"$1"-"$3"/("$7","$8","$9","$10")"}' ibdcheck02.rel | sort -k1,1 > relationships.rel

sqlite3 << EOF
.output relationships_merged.csv
.mode csv 
create table data(a,b);
.separator ' '
.import relationships.rel data
select a, group_concat(b,"/") from data group by a;
EOF


gzip -f ibdgenome.genome

echo -e '\n\n'
echo -e ${hru}
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'


echo '
\newpage
\section{IBD analysis of \texttt{\detokenize{'${dataset##*/}'}}}

In this section we describe the IBD analysis. The IBD analysis preceeds the Principal Component Analysis, and is used to detect high-IBD pairs. The cut-off for the pairs is a pi-hat of $0.1875$. We remove one sample of each pair, focusing on removing the least cases, and the samples that aggregate the most IBD pairings. For the IBD analysis, we first run an LD-pruning step.\\

\begin{table}[H]
	\caption{Overview of samples in the analysis}
	\centering
	\begin{tabular}{lllll}\toprule
		Dataset&Cases&Controls&Total&SNPs \\ \midrule
		Input dataset&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&
		'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&
		'`wc -l ${dataset}.fam | awk '{print $1}'`'&
		'`wc -l ${dataset}.bim | awk '{print $1}'`'\\
		IBD-ready dataset&'`awk '{if($6==2) print}' ibd.fam | wc -l`'&
		'`awk '{if($6==1) print}' ibd.fam | wc -l`'&
		'`wc -l ibd.fam | awk '{print $1}'`'&
		'`wc -l ibd.bim | awk '{print $1}'`' \\
		\bottomrule
	\end{tabular}
\end{table}

We show a histogram of the pair frequencies of all the IBD rates.\\

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e hist-plot-ibdgenome.genome_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the pairs generated by the IBD analysis}
\end{figure}

The next plots show the Z0,Z1 and Z2 scatter plots. The color represents the pihat $(=Z1+1/2*Z2)$

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e ibd-plot-ibdcheck01-6:7-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e ibd-plot-ibdcheck01-7:8-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}

\lstinputlisting[breaklines, 
				firstline=1, 
				frame=single, 
				title={Sample summary in the dataset after Individual QC}
						]{'`readlink -e ibd_check_summary.txt`'}

' >> report.tex

echo '
\section{IBD analysis of \texttt{\detokenize{'${dataset##*/}'}}}

\begin{frame}{IBD check of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e hist-plot-ibdgenome.genome_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the pairs generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e ibd-plot-ibdcheck01-6:7-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`\
		readlink -e ibd-plot-ibdcheck01-7:8-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \texttt{\detokenize{'${dataset##*/}'}}}
\lstinputlisting[breaklines, 
				firstline=1, 
				frame=single, 
				title={Sample summary in the dataset after Individual QC}
						]{'`readlink -e ibd_check_summary.txt`'}
\end{frame}
' >> slides.tex

echo '\end{document}' >> report.tex
pdflatex -interaction=batchmode report.tex
rm -rf report.aux report.log

echo '\end{document}' >> slides.tex
pdflatex -interaction=batchmode slides.tex
rm -rf slides.aux slides.log slides.nav slides.out slides.snm slides.toc
echo -e ${hrb}'\n'

echo -e $endheader

date >> runtime.info

cd ..
