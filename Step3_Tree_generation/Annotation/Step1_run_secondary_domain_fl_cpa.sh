#!/bin/bash
#SBATCH --job-name=otherdomains
#SBATCH --time=72:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=10GB          # Memory in MB
#SBATCH --cpus-per-task=10
#SBATCH --array=0-2
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/Other_domainRun_python_output%A-%a.out
#SBATCH --error=slurm_output/Other_domainRun_python_error%A-%a.err

#input are the fl sequences of the cpa homologs. 

declare -a array=("Eukarya" "Archaea" "Bacteria")
module load Python/3.11.6-foss-2023a

module load HMMER/3.3.2-GCC-12.3.0

IN=/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set
HMM_DB=/nesi/nobackup/uc04105/results/HMM
#NhaA
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF06965_NhaA.sthk ${HMM_DB}/PF06965.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF06965_NhaA.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF06965_NhaA 262

#cNMP
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF00027_cNMP.sthk ${HMM_DB}/PF00027.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF00027_cNMP.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF00027_cNMP 62

#CHX2nd
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23256_CHX2nd.sthk ${HMM_DB}/PF23256.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23256_CHX2nd.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23256_CHX2nd 99

#CHXC
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23259_CHXC.sthk ${HMM_DB}/PF23259.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23259_CHXC.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF23259_CHXC 105

#NhaC1term
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF08619_NhaC1term.sthk ${HMM_DB}/PF08619.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF08619_NhaC1term.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF08619_NhaC1term 325

#TrKAC
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02080_TrkC.sthk ${HMM_DB}/PF02080.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02080_TrkC.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02080_TrkC 50

#TrkAN
hmmalign -o ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02254_TrkN.sthk ${HMM_DB}/PF02254.hmm ${IN}/${array[$SLURM_ARRAY_TASK_ID]}_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta
python parse_stockholm_special.py ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02254_TrkN.sthk ${IN}/HMMalign/${array[$SLURM_ARRAY_TASK_ID]}_fl_vsPF02254_TrkN 82
