(use-modules (ice-9 threads))
(use-modules (srfi srfi-1))
;; Crisp rules about translating a link into another link

;; Helpers
(define ConceptT (TypeInh "ConceptNode"))
(define GeneT (Type "GeneNode"))

(define-public (gen-present-link-translation-rule LINK-1 LINK-TYPE-2 VAR-TYPE)
  (let* ((X (gar LINK-1))
         (Y (gdr LINK-1))
         (XY-2 (LINK-TYPE-2 X Y)))
    (Bind
      (VariableList
        (TypedVariable X VAR-TYPE)
        (TypedVariable Y VAR-TYPE))
      (Present
        XY-1)
      XY-2)))

(define-public (inheritance->subset)
    ;; get patient atoms and run the deduction in batch
    (cog-logger-info "Inheritance->Subset")
    ;;apply fc to get the relationship between go's and patients
    (let* ((atoms (cog-get-atoms 'InheritanceLink))
            (batch-num 0)
            (batch-size (/ (len genes) (current-processor-count)))
            (batch-ls (split-lst genes batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (ln)
                  (gen-present-link-translation-rule ln SubsetLink ConceptT)) batch)) batches)
        (cog-logger-info "Done!")))

(define-public (take-custom lst n)
    (if (< (length lst) n)
        (take lst (length lst))
        (take lst n)))

(define-public (drop-custom lst n)
    (if (< (length lst) n)
        (drop lst (length lst))
        (drop lst n)))
(define-public (split-lst lst n)
    (if (null? lst) '()
        (cons (take-custom lst n) (split-lst (drop-custom lst n) n))))
  
