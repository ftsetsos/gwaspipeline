#!/bin/sh
source ~/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2018 Jan 2
######---------------------------------------------------------



######---------------------------------------------------------
######    Start the log
######---------------------------------------------------------

date | tee -a ~/scripthistory_gwas.log
echo 'Running pca.sh with the following arguments:' \
		| tee -a ~/scripthistory_gwas.log

######---------------------------------------------------------

######---------------------------------------------------------
######    Parse arguments and write them on the log
######---------------------------------------------------------

OPTIND=1
dataset=''
outliers=${BLANKFILE}
projectfile='no'
threads=8
memory=16000

while getopts "hd:o:p:t:m:" opt; do
	case "$opt" in
		h) echo "Runs PCA on the data. Prunes it first."
			echo "Usage: pca.sh -d [dataset] "
			echo "			-o [outliers file (optional)] "
			echo "			-p [projection file (optional)] "
			echo "			-t [threads (optional: default 8)] "
			echo "			-m [memory (optional: default 16000)]"
			sed -i '$d' ~/scripthistory_gwas.log
			sed -i '$d' ~/scripthistory_gwas.log
			exit 0 ;;
		d) dataset=`readlink -e $OPTARG.bed | sed 's/\.bed//'` 
			echo -e '\t -d '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		o) outliers=`readlink -e $OPTARG` 
			echo -e '\t -o '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
		p) projectfile=`readlink -e $OPTARG` 
			echo -e '\t -p '$OPTARG | tee -a ~/scripthistory_gwas.log ;;
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

mkdir ${dataset##*/}_pca_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_pca_$(echo ${outliers##*/} | cut -d'.' -f1)

date >> runtime.info

echo -e $latexdoc > report.tex
echo -e $latexbeamer > slides.tex

echo -e $startheader




commandvar="awk '{if((\$1==$mhcchrom && \$4>$mhcstart && \$4<$mhcend) || \
					(\$1==$chr8invchrom && \$4>$chr8invstart && \$4<$chr8invstop) || \
					(\$1==$chr17invchrom && \$4>$chr17invstart && \$4<$chr17invstop)) \
						print \$2, \$1, \$4}' ${dataset}.bim > mhc817.snps"
eval $commandvar
awk '{if($1>22 || $1==0 ) print}' ${dataset}.bim > nonautosomal.snps
awk '{print $2}' ${dataset}.bim | awk '{print length, $0}' | awk '{if($1>39) print $2}' \
	> longname.snps
cat mhc817.snps nonautosomal.snps longname.snps > mhc817nonauto.snps


if [ "${projectfile}" != "no" ] ; then
	awk '{print $1,$2}' ${projectfile} > project.ind
	awk '{print $1,$2,3}' ${projectfile} > project.pheno
	cp ${dataset}.fam projected.pheno1
	while read line ; do
		sed -i "/$line/d" projected.pheno1
	done < project.ind
	awk '{print $1,$2,$6}' projected.pheno1 > projected.pheno
	cat project.pheno projected.pheno > all.pheno
else
	awk '{print $1,$2,$6}' ${dataset}.fam > all.pheno
fi

echo -e '\nStep ❶'
echo -e ${hru}
	${PLINK2PATH} --allow-no-sex \
					--bfile ${dataset} \
					--geno ${snpmissing} \
					--maf ${maf} \
					--hwe ${hwe} \
					--exclude mhc817nonauto.snps \
					--indep-pairwise 200 100 0.2 \
					--remove ${outliers} \
					--out pca1 \
					--memory ${memory}
echo -e ${hrb}'\n'

############   If there are still more than 150K SNPs in the dataset,
############   then prune again using the same settings.
if [ "$(wc -l pca1.prune.in | awk '{print $1}')" -gt 150000 ] ; then 

	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --allow-no-sex \
						--bfile ${dataset} \
						--extract pca1.prune.in \
						--indep-pairwise 200 100 0.2 \
						--remove ${outliers} \
						--out pca2 \
						--memory ${memory}
	echo -e ${hrb}'\n'
	
	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --allow-no-sex \
						--bfile ${dataset} \
						--pheno all.pheno \
						--extract pca2.prune.in \
						--make-bed \
						--remove ${outliers} \
						--out pca \
						--memory ${memory}
	echo -e ${hrb}'\n'
else
	echo -e '\nStep ❷'
	echo -e ${hru}
		${PLINK2PATH} --allow-no-sex \
						--bfile ${dataset} \
						--pheno all.pheno \
						--extract pca1.prune.in \
						--make-bed \
						--out pca \
						--remove ${outliers} \
						--memory ${memory}
	echo -e ${hrb}'\n'
fi

echo 'genotypename: pca.bed
snpname:      pca.bim
indivname:    pca.fam
evecoutname:  pca.evec
evaloutname:  pca.eval
numoutevec:   20
numoutlieriter: 0
numthreads: '${threads} > pca.parfile
if [ "$projectfile" != "no" ] ; then
echo '3' > projectpop.txt
echo 'poplistname: projectpop.txt' >> pca.parfile
fi
echo 'familynames: NO
altnormstyle: NO' >> pca.parfile


echo -e '\nStep ❸'
echo -e ${hru}
	${EIGPATH}bin/smartpca -p pca.parfile | tee pca.pcalog
echo -e ${hrb}''


echo -e '\n\n'
echo -e ${hru}
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 1:2
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 2:3
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 3:4
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 4:5
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 5:6
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 7:8
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 9:10
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 1:2:3
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 3:4:5
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 6:7:8
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 1:2:3
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 4:5:6
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 7:8:9
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20
echo -e ${hrb}''

cat <(printf 'FID IID ' && printf 'PC%s ' {1..20} && echo 'Status') \
	<(tail -n +2 pca.evec | sort -k1,1 | awk '{print $1,$0}') \
																> pca.covar

cat <(printf 'FID IID ' && printf 'PC%s ' {1..20} && echo 'Status Sex' ) \
	<(join -1 1 -2 1 <(awk '{print $2,$5}' pca.fam | sort -k1,1 ) \
					<(sort -k1,1 pca.covar) \
					| awk '{print $0,$2}' | cut -d ' ' -f1,3- ) > pcasex.covar

echo -e ${hru}
echo 'Running MatchIt on the pca data.'
Rscript ${SCRIPTPATH}samplematching.R
echo -e ${hrb}

echo -e '\n\n'
echo -e ${hru}
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'


echo '
\newpage
\section{PCA of \texttt{\detokenize{'${dataset##*/}'}}}

In this section we describe the results of the Principal Component Analysis. The PCA step is run after the IBD analysis, removing the IBD individuals that could lead to confounding. The PCA step is also rerun a second time to remove population outliers from the PCA. For the PCA, we first run an LD-pruning step.

\begin{table}[H]
	\caption{Overview of samples in the analysis}
	\centering
	\begin{tabular}{lllll}\toprule
		Dataset&Cases&Controls&Total&SNPs \\ \midrule
		Input dataset&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&
		'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&
		'`wc -l ${dataset}.fam | awk '{print $1}'`'&
		'`wc -l ${dataset}.bim | awk '{print $1}'`'\\
		PCA-ready dataset&'`awk '{if($6==2) print}' pca.fam | wc -l`'&
		'`awk '{if($6==1) print}' pca.fam | wc -l`'&
		'`wc -l pca.fam | awk '{print $1}'`'&
		'`wc -l pca.bim | awk '{print $1}'`' \\
		\bottomrule
	\end{tabular}
\end{table}

\begin{table}[H]
	\caption{Eigenvalues of the eigenvectors generated}
	\centering
	\begin{tabular}{l|llllllllll}\toprule
		Eigenvectors&1&2&3&4&5&6&7&8&9&10 \\ \midrule
		Eigenvalues&'`head -1 pca.evec | xargs | cut -d" " -f2-11 | tr " " "&" `' \\
		\bottomrule
	\end{tabular}
\end{table}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-1:2-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first two PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-3:4-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the third and fourth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-5:6-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the fifth and sixth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-7:8-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the seventh and eighth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-9:10-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the ninth and tenth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-1:2:3-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-3:4:5-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-1:2:3-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-4:5:6-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-7:8:9-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}
' >> report.tex

echo '
\section{PCA of \texttt{\detokenize{'${dataset##*/}'}}}

\begin{frame}{PCA of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-1:2-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first two PCs generated by the PCA}
\end{figure}
\end{frame}

\begin{frame}{PCA of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-3:4-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the third and fourth PCs generated by the PCA}
\end{figure}
\end{frame}

\begin{frame}{PCA of \texttt{\detokenize{'${dataset##*/}'}}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`\
		readlink -e pca-plot-pca-1:2:3-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}
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
