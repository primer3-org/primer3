;; Clisp binding for libprimer3
;; kumi (kkhawi@gmail.com) Dec-16-2006  
;; Please don't overwrite this file with a newly generated one. Swig support for clisp is not in good shape
;; and there are a few bugs. This file is *heavily* edited and hand tweaked for hours.
;; 
(defpackage :libprimer3
  (:use :common-lisp :ffi)
  (:export
	:make-program_args
	:program_args-format_output
	:program_args-twox_compat
	:program_args-strict_tags
	:make-seq_lib
	:seq_lib-names
	:seq_lib-seqs
	:seq_lib-rev_compl_seqs
	:seq_lib-weight
	:seq_lib-repeat_file
	:seq_lib-error
	:seq_lib-warning
	:seq_lib-seq_num
	:make-oligo_weights
	:oligo_weights-temp_gt
	:oligo_weights-temp_lt
	:oligo_weights-gc_content_gt
	:oligo_weights-gc_content_lt
	:oligo_weights-compl_any
	:oligo_weights-compl_end
	:oligo_weights-repeat_sim
	:oligo_weights-length_lt
	:oligo_weights-length_gt
	:oligo_weights-seq_quality
	:oligo_weights-end_quality
	:oligo_weights-pos_penalty
	:oligo_weights-end_stability
	:oligo_weights-num_ns
	:make-pair_weights
	:pair_weights-primer_quality
	:pair_weights-io_quality
	:pair_weights-diff_tm
	:pair_weights-compl_any
	:pair_weights-compl_end
	:pair_weights-product_tm_lt
	:pair_weights-product_tm_gt
	:pair_weights-product_size_lt
	:pair_weights-product_size_gt
	:pair_weights-repeat_sim
	:make-primer_args
	:primer_args-pr_min
	:primer_args-pr_max
	:primer_args-repeat_lib
	:primer_args-io_mishyb_library
	:primer_args-primer_weights
	:primer_args-io_weights
	:primer_args-pr_pair_weights
	:primer_args-glob_err
	:primer_args-opt_tm
	:primer_args-min_tm
	:primer_args-max_tm
	:primer_args-max_diff_tm
	:primer_args-opt_gc_content
	:primer_args-max_gc
	:primer_args-min_gc
	:primer_args-salt_conc
	:primer_args-dna_conc
	:primer_args-io_opt_tm
	:primer_args-io_min_tm
	:primer_args-io_max_tm
	:primer_args-io_opt_gc_content
	:primer_args-io_max_gc
	:primer_args-io_min_gc
	:primer_args-io_salt_conc
	:primer_args-io_dna_conc
	:primer_args-outside_penalty
	:primer_args-inside_penalty
	:primer_args-product_max_tm
	:primer_args-product_min_tm
	:primer_args-product_opt_tm
	:primer_args-max_end_stability
	:primer_args-num_intervals
	:primer_args-num_ns_accepted
	:primer_args-primer_task
	:primer_args-file_flag
	:primer_args-explain_flag
	:primer_args-primer_opt_size
	:primer_args-primer_min_size
	:primer_args-primer_max_size
	:primer_args-product_opt_size
	:primer_args-io_num_ns_accepted
	:primer_args-io_primer_opt_size
	:primer_args-io_primer_min_size
	:primer_args-io_primer_max_size
	:primer_args-gc_clamp
	:primer_args-liberal_base
	:primer_args-max_poly_x
	:primer_args-io_max_poly_x
	:primer_args-first_base_index
	:primer_args-num_return
	:primer_args-min_quality
	:primer_args-min_end_quality
	:primer_args-quality_range_min
	:primer_args-quality_range_max
	:primer_args-io_min_quality
	:primer_args-io_min_end_quality
	:primer_args-pick_anyway
	:primer_args-repeat_compl
	:primer_args-io_repeat_compl
	:primer_args-pair_repeat_compl
	:primer_args-self_any
	:primer_args-self_end
	:primer_args-io_self_any
	:primer_args-io_self_end
	:primer_args-pair_compl_any
	:primer_args-pair_compl_end
	:make-pair_stats
	:pair_stats-considered
	:pair_stats-product
	:pair_stats-target
	:pair_stats-temp_diff
	:pair_stats-compl_any
	:pair_stats-compl_end
	:pair_stats-internal
	:pair_stats-repeat_sim
	:pair_stats-high_tm
	:pair_stats-low_tm
	:pair_stats-ok
	:make-oligo_stats
	:oligo_stats-considered
	:oligo_stats-ns
	:oligo_stats-target
	:oligo_stats-excluded
	:oligo_stats-gc
	:oligo_stats-gc_clamp
	:oligo_stats-temp_min
	:oligo_stats-temp_max
	:oligo_stats-compl_any
	:oligo_stats-compl_end
	:oligo_stats-repeat
	:oligo_stats-poly_x
	:oligo_stats-seq_quality
	:oligo_stats-stability
	:oligo_stats-no_orf
	:oligo_stats-ok
	:make-seq_args
	:seq_args-error
	:seq_args-warning
	:seq_args-num_targets
	:seq_args-tar
	:seq_args-num_excl
	:seq_args-excl
	:seq_args-num_internal_excl
	:seq_args-excl_internal
	:seq_args-incl_s
	:seq_args-incl_l
	:seq_args-start_codon_pos
	:seq_args-stop_codon_pos
	:seq_args-quality
	:seq_args-sequence
	:seq_args-sequence_name
	:seq_args-sequence_file
	:seq_args-trimmed_seq
	:seq_args-left_input
	:seq_args-right_input
	:seq_args-internal_input
	:seq_args-left_expl
	:seq_args-right_expl
	:seq_args-intl_expl
	:seq_args-pair_expl
	:make-dpal_args
	:dpal_args-check_chars
	:dpal_args-debug
	:dpal_args-fail_stop
	:dpal_args-flag
	:dpal_args-force_generic
	:dpal_args-force_long_generic
	:dpal_args-force_long_maxgap1
	:dpal_args-gap
	:dpal_args-gapl
	:dpal_args-max_gap
	:dpal_args-score_max
	:dpal_args-score_only
	:dpal_args-ssm
	:make-primer_pair
	:primer_pair-pair_quality
	:primer_pair-compl_measure
	:primer_pair-diff_tm
	:primer_pair-product_tm
	:primer_pair-product_tm_oligo_tm_diff
	:primer_pair-t_opt_a
	:primer_pair-compl_any
	:primer_pair-compl_end
	:primer_pair-repeat_sim
	:primer_pair-left
	:primer_pair-right
	:primer_pair-intl
	:primer_pair-product_size
	:primer_pair-target
	:primer_pair-rep_name
	:make-pair_array_t
	:pair_array_t-storage_size
	:pair_array_t-num_pairs
	:pair_array_t-pairs
	:make-primer_state
	:primer_state-local_args
	:primer_state-local_end_args
	:primer_state-end_args
	:primer_state-local_args_ambig
	:primer_state-local_end_args_ambig
	:primer_state-f
	:primer_state-r
	:primer_state-mid
	:primer_state-n_f
	:primer_state-n_r
	:primer_state-n_m
	:primer_state-f_len
	:primer_state-r_len
	:primer_state-mid_len
	:primer_state-best_pairs
	:primer_state-err
	:PR_ERR_NONE
	:PR_ERR_OUT_OF_MEMORY
	:PR_ERR_CANNOT_OPEN_FILE
	:PR_ERR_ALIGNMENT_FAILED
	:make-primer_error
	:primer_error-system_errno
	:primer_error-local_errno
	:primer_error-error_msg
	:primer_error-jmpenv
	:OV_UNINITIALIZED
	:OV_OK
	:OV_TOO_MANY_NS
	:OV_INTERSECT_TARGET
	:OV_GC_CONTENT
	:OV_TM_LOW
	:OV_TM_HIGH
	:OV_SELF_ANY
	:OV_SELF_END
	:OV_EXCL_REGION
	:OV_GC_CLAMP
	:OV_END_STAB
	:OV_POLY_X
	:OV_SEQ_QUALITY
	:OV_LIB_SIM
	:make-rep_sim
	:rep_sim-name
	:rep_sim-min
	:rep_sim-max
	:rep_sim-score
	:make-primer_rec
	:primer_rec-repeat_sim
	:primer_rec-temp
	:primer_rec-gc_content
	:primer_rec-position_penalty
	:primer_rec-quality
	:primer_rec-end_stability
	:primer_rec-start
	:primer_rec-seq_quality
	:primer_rec-self_any
	:primer_rec-self_end
	:primer_rec-target
	:primer_rec-excl
	:primer_rec-ok
	:primer_rec-length
	:primer_rec-num_ns
	:primer_rec-position_penalty_infinite
	:primer_rec-must_use
	:make-pr_append_str
	:pr_append_str-storage_size
	:pr_append_str-data
	:OT_LEFT
	:OT_RIGHT
	:OT_INTL))

(in-package :libprimer3)

(default-foreign-language :stdc)


(ffi:def-c-struct program_args
		  (format_output character)
		  (twox_compat character)
		  (strict_tags character))

(ffi:def-c-type program_args program_args)

(ffi:def-c-struct pr_append_str
		  (storage_size ffi:int)
		  (data ffi:c-string))

(ffi:def-c-type pr_append_str pr_append_str)

(ffi:def-c-struct seq_lib
		  (names (ffi:c-ptr ffi:c-string))
		  (seqs (ffi:c-ptr ffi:c-string))
		  (rev_compl_seqs (ffi:c-ptr ffi:c-string))
		  (weight (ffi:c-ptr DOUBLE-FLOAT))
		  (repeat_file ffi:c-string)
		  (error pr_append_str)
		  (warning pr_append_str)
		  (seq_num ffi:int))

(ffi:def-c-type seq_lib seq_lib)

(ffi:def-c-struct oligo_weights
		  (temp_gt DOUBLE-FLOAT)
		  (temp_lt DOUBLE-FLOAT)
		  (gc_content_gt DOUBLE-FLOAT)
		  (gc_content_lt DOUBLE-FLOAT)
		  (compl_any DOUBLE-FLOAT)
		  (compl_end DOUBLE-FLOAT)
		  (repeat_sim DOUBLE-FLOAT)
		  (length_lt DOUBLE-FLOAT)
		  (length_gt DOUBLE-FLOAT)
		  (seq_quality DOUBLE-FLOAT)
		  (end_quality DOUBLE-FLOAT)
		  (pos_penalty DOUBLE-FLOAT)
		  (end_stability DOUBLE-FLOAT)
		  (num_ns DOUBLE-FLOAT))

(ffi:def-c-type oligo_weights oligo_weights)


(defconstant PR_MAX_INTERVAL_ARRAY 200)

(ffi:def-c-struct oligo_weigts
		  (temp_gt DOUBLE-FLOAT)
		  (temp_lt DOUBLE-FLOAT)
		  (gc_content_gt DOUBLE-FLOAT)
		  (gc_content_lt DOUBLE-FLOAT)
		  (compl_any DOUBLE-FLOAT)
		  (compl_end DOUBLE-FLOAT)
		  (repeat_sim DOUBLE-FLOAT)
		  (length_lt DOUBLE-FLOAT)
		  (length_gt DOUBLE-FLOAT)
		  (seq_quality DOUBLE-FLOAT)
		  (end_quality DOUBLE-FLOAT)
		  (pos_penalty DOUBLE-FLOAT)
		  (end_stability DOUBLE-FLOAT)
		  (num_ns DOUBLE-FLOAT))

(ffi:def-c-type oligo_weigts oligo_weigts)

(ffi:def-c-struct pair_weights
		  (primer_quality DOUBLE-FLOAT)
		  (io_quality DOUBLE-FLOAT)
		  (diff_tm DOUBLE-FLOAT)
		  (compl_any DOUBLE-FLOAT)
		  (compl_end DOUBLE-FLOAT)
		  (product_tm_lt DOUBLE-FLOAT)
		  (product_tm_gt DOUBLE-FLOAT)
		  (product_size_lt DOUBLE-FLOAT)
		  (product_size_gt DOUBLE-FLOAT)
		  (repeat_sim DOUBLE-FLOAT))

(ffi:def-c-type pair_weigts pair_weights)

(ffi:def-c-enum task 
		(pick_pcr_primers 0)
		(pick_pcr_primers_and_hyb_probe 1)
		(pick_left_only 2)(pick_right_only 3)
		(pick_hyb_probe_only 4))

(ffi:def-c-type task task)


(ffi:def-c-struct primer_args
		  (pr_min (ffi:c-array ffi:int 200))
		  (pr_max (ffi:c-array ffi:int 200))
		  (repeat_lib seq_lib)
		  (io_mishyb_library seq_lib)
		  (primer_weights oligo_weights)
		  (io_weights oligo_weights)
		  (pr_pair_weights pair_weights)
		  (glob_err pr_append_str)
		  (opt_tm DOUBLE-FLOAT)
		  (min_tm DOUBLE-FLOAT)
		  (max_tm DOUBLE-FLOAT)
		  (max_diff_tm DOUBLE-FLOAT)
		  (opt_gc_content DOUBLE-FLOAT)
		  (max_gc DOUBLE-FLOAT)
		  (min_gc DOUBLE-FLOAT)
		  (salt_conc DOUBLE-FLOAT)
		  (dna_conc DOUBLE-FLOAT)
		  (io_opt_tm DOUBLE-FLOAT)
		  (io_min_tm DOUBLE-FLOAT)
		  (io_max_tm DOUBLE-FLOAT)
		  (io_opt_gc_content DOUBLE-FLOAT)
		  (io_max_gc DOUBLE-FLOAT)
		  (io_min_gc DOUBLE-FLOAT)
		  (io_salt_conc DOUBLE-FLOAT)
		  (io_dna_conc DOUBLE-FLOAT)
		  (outside_penalty DOUBLE-FLOAT)
		  (inside_penalty DOUBLE-FLOAT)
		  (product_max_tm DOUBLE-FLOAT)
		  (product_min_tm DOUBLE-FLOAT)
		  (product_opt_tm DOUBLE-FLOAT)
		  (max_end_stability DOUBLE-FLOAT)
		  (num_intervals ffi:int)
		  (num_ns_accepted ffi:int)
		  (primer_task task)
		  (file_flag ffi:int)
		  (explain_flag ffi:int)
		  (primer_opt_size ffi:int)
		  (primer_min_size ffi:int)
		  (primer_max_size ffi:int)
		  (product_opt_size ffi:int)
		  (io_num_ns_accepted ffi:int)
		  (io_primer_opt_size ffi:int)
		  (io_primer_min_size ffi:int)
		  (io_primer_max_size ffi:int)
		  (gc_clamp ffi:int)
		  (liberal_base ffi:int)
		  (max_poly_x ffi:int)
		  (io_max_poly_x ffi:int)
		  (first_base_index ffi:int)
		  (num_return ffi:int)
		  (min_quality ffi:int)
		  (min_end_quality ffi:int)
		  (quality_range_min ffi:int)
		  (quality_range_max ffi:int)
		  (io_min_quality ffi:int)
		  (io_min_end_quality ffi:int)
		  (pick_anyway ffi:int)
		  (repeat_compl ffi:short)
		  (io_repeat_compl ffi:short)
		  (pair_repeat_compl ffi:short)
		  (self_any ffi:short)
		  (self_end ffi:short)
		  (io_self_any ffi:short)
		  (io_self_end ffi:short)
		  (pair_compl_any ffi:short)
		  (pair_compl_end ffi:short))

(ffi:def-c-type primerargs primer_args)

(ffi:def-c-type interval_array_t (ffi:c-array ffi:int (200 2))) ;PR_MAX_INTERVAL_ARRAY

(ffi:def-c-type dpal_ssm (ffi:c-array ffi:int (256 256))) ; each UCHAR_MAX + 1

(ffi:def-c-struct pair_stats
		  (considered ffi:int)
		  (product ffi:int)
		  (target ffi:int)
		  (temp_diff ffi:int)
		  (compl_any ffi:int)
		  (compl_end ffi:int)
		  (internal ffi:int)
		  (repeat_sim ffi:int)
		  (high_tm ffi:int)
		  (low_tm ffi:int)
		  (ok ffi:int))

(ffi:def-c-type pair_stats pair_stats)

(ffi:def-c-struct oligo_stats
		  (considered ffi:int)
		  (ns ffi:int)
		  (target ffi:int)
		  (excluded ffi:int)
		  (gc ffi:int)
		  (gc_clamp ffi:int)
		  (temp_min ffi:int)
		  (temp_max ffi:int)
		  (compl_any ffi:int)
		  (compl_end ffi:int)
		  (repeat ffi:int)
		  (poly_x ffi:int)
		  (seq_quality ffi:int)
		  (stability ffi:int)
		  (no_orf ffi:int)
		  (ok ffi:int))

(ffi:def-c-type oligo_stats oligo_stats)

(ffi:def-c-struct seq_args
		  (error pr_append_str)
		  (warning pr_append_str)
		  (num_targets ffi:int)
		  (tar interval_array_t)
		  (num_excl ffi:int)
		  (excl interval_array_t)
		  (num_internal_excl ffi:int)
		  (excl_internal interval_array_t)
		  (incl_s ffi:int)
		  (incl_l ffi:int)
		  (start_codon_pos ffi:int)
		  (stop_codon_pos ffi:int)
		  (quality (ffi:c-ptr ffi:int))
		  (sequence ffi:c-string)
		  (sequence_name ffi:c-string)
		  (sequence_file ffi:c-string)
		  (trimmed_seq ffi:c-string)
		  (left_input ffi:c-string)
		  (right_input ffi:c-string)
		  (internal_input ffi:c-string)
		  (left_expl oligo_stats)
		  (right_expl oligo_stats)
		  (intl_expl oligo_stats)
		  (pair_expl pair_stats))


(ffi:def-c-type seqargs seq_args)

(ffi:def-c-struct dpal_args
		  (check_chars ffi:int)
		  (debug ffi:int)
		  (fail_stop ffi:int)
		  (flag ffi:int)
		  (force_generic ffi:int)
		  (force_long_generic ffi:int)
		  (force_long_maxgap1 ffi:int)
		  (gap ffi:int)
		  (gapl ffi:int)
		  (max_gap ffi:int)
		  (score_max ffi:int)
		  (score_only ffi:int)
		  (ssm dpal_ssm))

(ffi:def-c-struct rep_sim
		  (name character)
		  (min ffi:short)
		  (max ffi:short)
		  (score (ffi:c-ptr ffi:short)))

(ffi:def-c-type rep_sim rep_sim)

(ffi:def-c-enum oligo_violation 
		(OV_UNINITIALIZED -1)
		(OV_OK 0)(OV_TOO_MANY_NS 1)
		(OV_INTERSECT_TARGET 2)(OV_GC_CONTENT 3)
		(OV_TM_LOW 4)(OV_TM_HIGH 5)(OV_SELF_ANY 6)
		(OV_SELF_END 7)(OV_EXCL_REGION 8)(OV_GC_CLAMP 9)
		(OV_END_STAB 10)(OV_POLY_X 11)(OV_SEQ_QUALITY 12)(OV_LIB_SIM 13))

(ffi:def-c-struct primer_rec
		  (repeat_sim rep_sim)
		  (temp DOUBLE-FLOAT)
		  (gc_content DOUBLE-FLOAT)
		  (position_penalty DOUBLE-FLOAT)
		  (quality DOUBLE-FLOAT)
		  (end_stability DOUBLE-FLOAT)
		  (start ffi:int)
		  (seq_quality ffi:int)
		  (self_any ffi:short)
		  (self_end ffi:short)
		  (target character)
		  (excl character)
		  (ok oligo_violation)
		  (length character)
		  (num_ns character)
		  (position_penalty_infinite character)
		  (must_use character))

(ffi:def-c-struct primer_pair
		  (pair_quality DOUBLE-FLOAT)
		  (compl_measure DOUBLE-FLOAT)
		  (diff_tm DOUBLE-FLOAT)
		  (product_tm DOUBLE-FLOAT)
		  (product_tm_oligo_tm_diff DOUBLE-FLOAT)
		  (t_opt_a DOUBLE-FLOAT)
		  (compl_any ffi:int)
		  (compl_end ffi:int)
		  (repeat_sim ffi:short)
		  (left (ffi:c-pointer primer_rec))
		  (right (ffi:c-pointer primer_rec))
		  (intl (ffi:c-pointer primer_rec))
		  (product_size ffi:int)
		  (target ffi:int)
		  (rep_name ffi:c-string))

(ffi:def-c-type primpair primer_pair)

(ffi:def-c-struct pair_array_t
		  (storage_size ffi:int)
		  (num_pairs ffi:int)
		  (pairs (ffi:c-pointer primer_pair)))

(ffi:def-c-type pair_array_t pair_array_t)

(ffi:def-c-enum primer_errno 
		(PR_ERR_NONE 0)(PR_ERR_OUT_OF_MEMORY 1)
		(PR_ERR_CANNOT_OPEN_FILE 2)(PR_ERR_ALIGNMENT_FAILED 3))

(ffi:def-c-type primer_errno primer_errno)

(ffi:def-c-struct primer_error
		  (system_errno ffi:int)
		  (local_errno primer_errno)
		  (error_msg ffi:c-string)
		  (jmpenv (ffi:c-ptr NIL)))

(ffi:def-c-type primer_error primer_error)

(ffi:def-c-struct primer_state
		  (local_args dpal_args)
		  (local_end_args dpal_args)
		  (end_args dpal_args)
		  (local_args_ambig dpal_args)
		  (local_end_args_ambig dpal_args)
		  (f (ffi:c-pointer primer_rec))
		  (r (ffi:c-pointer primer_rec))
		  (mid (ffi:c-pointer primer_rec))
		  (n_f ffi:int)
		  (n_r ffi:int)
		  (n_m ffi:int)
		  (f_len ffi:int)
		  (r_len ffi:int)
		  (mid_len ffi:int)
		  (best_pairs pair_array_t)
		  (err primer_error))

(ffi:def-c-type primer_state primer_state)


; who uses this?
(ffi:def-c-enum oligo_type (OT_LEFT 0)(OT_RIGHT 1)(OT_INTL 2))

(ffi:def-c-type oligo_type oligo_type)

