#!/bin/sh
source ~/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2018 Jan 2
######---------------------------------------------------------



######---------------------------------------------------------
######    Start the log
######---------------------------------------------------------

date | tee -a ~/scripthistory_gwas.log
echo 'Running prsice.sh with the following arguments:' \
		| tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------

######---------------------------------------------------------
######    Parse arguments and write them on the log
######---------------------------------------------------------

OPTIND=1
dataset=''
threads=12
outliers=${BLANKFILE}
format='plink'
covarfile='none'
covarcol='@PC[1-5],Sex'
prevalence='0.01'
memory=50000

while getopts "hd:o:s:c:n:f:t:" opt; do
	case "$opt" in
		h) echo "Performs polygenic risk scoring using PRSice."
			echo "Usage: prsice.sh -d [dataset] "
			echo "			-s [summary statistics file]"
			echo "			-o [individuals to remove (default: null)]"
			echo "			-c [covariate file (default: none)]"
			echo "			-k [prevalence of the disorder (default: 0.01)]"
			echo "			-f [summary statistics file format (default: plink)]"
			echo "			-t [threads (optional: default 1)] "
			sed -i '$d' ~/scripthistory_gwas.log
			sed -i '$d' ~/scripthistory_gwas.log
			exit 0 ;;
		d) dataset=`readlink -e $OPTARG.bed | sed 's/\.bed//'` 
			echo -e '\t -d '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		o) outliers=`readlink -e $OPTARG` 
			echo -e '\t -o '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		s) sumstats=`readlink -e $OPTARG` 
			echo -e '\t -s '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		c) covarfile=`readlink -e $OPTARG` 
			echo -e '\t -c '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		k) prevalence=`readlink -e $OPTARG` 
			echo -e '\t -k '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		n) covarcol=$OPTARG 
			echo -e '\t -n '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		f) format=$OPTARG 
			echo -e '\t -f '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		t) threads=$OPTARG 
			echo -e '\t -t '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
	esac
done
echo '' | tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------


######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_prs_${sumstats##*/}_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_prs_${sumstats##*/}_$(echo ${outliers##*/} | cut -d'.' -f1)

echo -e $latexdoc > report.tex
echo -e $latexbeamer > slides.tex

echo -e $startheader



#######---------------------------------------------------------
#######      First remove outliers or duplicates from the list
#######---------------------------------------------------------

echo -e '\nStep ❶'
echo -e ${hru}
${PLINK2PATH} --bfile ${dataset} \
				--allow-no-sex \
				--remove ${outliers} \
				--make-bed \
				--out subset \
				--threads ${threads} \
				--memory ${memory}
echo -e ${hrb}'\n'


#######---------------------------------------------------------
#######      Perform PRS based on the format of the sumstats
#######---------------------------------------------------------

if [ "$covarfile" != 'none' ]; then 
	covararg='-C '${covarfile}' -c '${covarcol}
fi

if [ "$format" == "plink" ] ; then
	echo -e '\nStep ❷'
	echo -e ${hru}
	Rscript ${PRSICEPATH}PRSice.R -b ${sumstats} \
								--chr CHR \
								--A1 A1 \
								--A2 A2 \
								--stat OR \
								--pvalue P \
								--snp SNP \
								--bp BP \
								--prevalence ${prevalence} \
								--target subset \
								--quantile 100 \
								--quant-break 1,5,10,20,40,60,80,90,95,99,100 \
								--quant-ref 60 \
								--binary-target T \
	 							$covararg \
								--perm 10000 \
								--thread ${threads} \
								--prsice ${PRSICEPATH}PRSice_linux
	echo -e ${hrb}'\n'
fi

join -1 1 -2 1 <(sort -k1,1 subset.fam) \
				<(sort -k1,1 PRSice.best) | awk '{print $6,$9}' > PRSice.groups

python3 ~/Dropbox/Fotis/myscripts/prsice_hist.py PRSice.groups

rm -rf subset.bed subset.bim subset.fam


#######---------------------------------------------------------
#######      Make LaTeX report
#######---------------------------------------------------------

echo '
\section{PRSice}

\begin{frame}{PRSice run data}
	\noindent Analysis ran on PRSice v'`head -1 PRSice.log | awk '{print $2}'`'\newline
	~\newline
	\noindent Number of variants included in the analysis: 
	'`grep 'total variant' PRSice.log | awk '{print $1}'`'\newline
	\noindent Clump window: '`grep clump-kb PRSice.log | awk '{print $2}'`'\newline
	\noindent Clump $r^2$: '`grep clump-r2 PRSice.log | awk '{print $2}'`'\newline
	\noindent Number of variants after clumping:
	'`grep 'clumping' PRSice.log | awk '{print $7}'`'\newline
	
	\begin{table}[H]
		\caption{Best p-value threshold calculated by PRS}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{lllllll}\toprule
			'`awk '{print $3,$4,$5,$6,$10,$11,$12}' PRSice.summary \
					| sed 's/Num_SNP/NumSNP/g' | sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
\end{frame}

\begin{frame}{PRSice}
\begin{figure}[H]
	\centering
	\includegraphics[height=0.7\textheight,keepaspectratio]{{'`echo PRSice_GROUPS* \
															| sed 's/\.png//' `'}.png}
	\caption{PRSice histogram}
\end{figure}
\end{frame}

\begin{frame}{PRSice}
\begin{figure}[H]
	\centering
	\includegraphics[height=0.8\textheight,keepaspectratio]{{'`echo PRSice_BARPLOT_* \
															| sed 's/\.png//' `'}.png}
	\caption{PRSice barplot}
\end{figure}
\end{frame}

\begin{frame}{PRSice}
\begin{figure}[H]
	\centering
	\includegraphics[height=0.8\textheight,keepaspectratio]{{'`echo PRSice_HIGH-RES_PLOT_* \
															| sed 's/\.png//' `'}.png}
	\caption{PRSice high-resolution p-value bins}
\end{figure}
\end{frame}

\begin{frame}{PRSice}
\begin{figure}[H]
	\centering
	\includegraphics[height=0.8\textheight,keepaspectratio]{{'`echo PRSice_STRATA_* \
															| sed 's/\.png//' `'}.png}
	\caption{PRSice quantile plot}
\end{figure}
\end{frame}
' >> slides.tex

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

echo -e $endheader

cd ..
