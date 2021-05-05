(define-module (pln-bio rules transitivity)
   #:use-module (opencog) 
   #:use-module (opencog exec) 
   #:use-module (opencog ure) 
   #:use-module (opencog logger)
   #:use-module (srfi srfi-1)
   #:use-module (pln-bio rules rule-utils)
)
;; Crisp rules about transitivity of predicates or inheritance links

;; Helpers
(define ConceptT (TypeInh "ConceptNode"))
(define GeneT (Type "GeneNode"))

(define-public (gen-present-link-transitivity-rule LINK LINK-TYPE VAR-TYPE)
  (let* ((X (gar LINK))
         (Y (gdr LINK))
         (Z (Variable "$Z"))
         (YZ (LINK-TYPE Y Z))
         (XZ (LINK-TYPE X Z)))
    (Bind
      (VariableList
        (TypedVariable X VAR-TYPE)
        (TypedVariable Y VAR-TYPE)
        (TypedVariable Z VAR-TYPE))
      (And
        (Present
          LINK
          YZ)
        (Not (Identical X Z)))
      XZ)))

(define-public (gen-present-mixed-link-transitivity-rule LINK-1 LINK-TYPE-1 LINK-TYPE-2
                                                  X-TYPE Y-TYPE Z-TYPE)
  (let* ((X (gar LINK-1))
         (Y (gdr LINK-1))
         (Z (Variable "$Z"))
         (YZ (LINK-TYPE-2 Y Z))
         (XZ (LINK-TYPE-1 X Z)))
    (Bind
      (VariableList
        (TypedVariable X X-TYPE)
        (TypedVariable Y Y-TYPE)
        (TypedVariable Z Z-TYPE))
      (And
        (Present
          LINK-1
          YZ)
        (Not (Identical X Z)))
      XZ)))

(define-public (gen-present-link-transitivity TYPE)
    (cog-logger-info "Running gen-present-link-transitivity")
    (let* ((atoms (cog-get-atoms TYPE))
            (batch-num 0)
            (batch-size (round (/ (length atoms) (current-processor-count))))
            (batch-ls (split-lst atoms batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (ln)
                  (gen-present-link-transitivity-rule ln TYPE ConceptT)) batch)) batches)
        (cog-logger-info "Done!")))

(define-public (gen-present-mixed-link-transitivity LINK-1-TYPE LINK-2-TYPE)
    ;; get patient atoms and run the deduction in batch
    (cog-logger-info "Running gen-present-link-transitivity")
    ;;apply fc to get the relationship between go's and patients
    (let* ((atoms (cog-get-atoms LINK-1-TYPE))
            (batch-num 0)
            (batch-size (round (/ (length atoms) (current-processor-count))))
            (batch-ls (split-lst atoms batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (ln)
                  (gen-present-mixed-link-transitivity-rule ln LINK-1-TYPE LINK-2-TYPE
                    GeneT ConceptT ConceptT)) batch)) batches)
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