(use-modules (pln-bio expr) (opencog logger) (pln-bio bio-utils)) 
;; install https://github.com/Habush/pln-brca-xp.git for these modules

(define (run-expr-deduction overexpr? path)
    (cog-logger-set-stdout! #t)
    (if overexpr?
        (begin 
            (cog-logger-info "Loading patient overexpression data")
            ;;load the atomese form of overexpr & underexpr
            (primitive-load (string-append path "/patient_gene_over_expr.scm"))      
            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (overexpression-dist)
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (get-overexpr-eval-ln)
            (write-atoms-to-file (string-append path "/normalized_patient_gene_over_expr.scm") (cog-get-atoms 'EvaluationLink))
            (cog-logger-info "Done!")
        )
        (begin 
            (cog-logger-info "Loading patient underexpression data")
            ;;load the atomese form of overexpr & underexpr
            (primitive-load (string-append path "/patient_gene_under_expr.scm"))
            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (underexpression-dist)
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (get-underexpr-eval-ln)
            (write-atoms-to-file (string-append path "/normalized_patient_gene_under_expr.scm") (cog-get-atoms 'EvaluationLink))
            (cog-logger-info "Done!"))))

(define path "/home/hedra/SIngularityNEt/knowledge-import/cancer/moses_ft83")
(run-expr-deduction #t path)
(clear)
(run-expr-deduction #f path)