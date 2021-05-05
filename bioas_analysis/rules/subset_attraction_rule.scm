(use-modules (ice-9 threads))
(use-modules (srfi srfi-1))

;; =============================================================================
;; Attraction introduction rule 
;;
;; Subset <STV>
;;   Set A
;;   B
;; Subset <SNTV>
;;   Not
;;     Set A
;;   B
;; |-
;; Attraction <TV>
;;   Set A
;;   B
;;
;; where TV is defined as follows
;;
;; TV.s = pattern-of(B, A)
;; TV.c = min(STV.c, SNTV.c)
;;
;; pattern-of(B,A) = (P(B|A)-P(B|¬A))+
;;
;; where s(B) is the prior of B and x+ is the positive part of x. For
;; now the prior of B is 1.

(define (subset-attraction-rule A)
  (let* ((B (Variable "$B")))
    (BindLink
      (Present
        (Subset (Set A) B)
        (Subset (Not (Set A)) B))
      (ExecutionOutputLink
        (GroundedSchemaNode "scm: gen-attraction-introduction")
        (ListLink
          ;; Conclusion
          (Attraction (Set A) B)
          ;; Premises
          (Subset (Set A) B)
          (Subset (Not (Set A)) B))))))

;; Formula
(define (gen-attraction-introduction conclusion . premises)
  (if (= (length premises) 2)
      (let* ((ATT conclusion)
             (SAB (car premises))
             (SNAB (cadr premises))
             (ATTs (max 0 (- (cog-mean SAB) (cog-mean SNAB))))
             (ATTc (min (cog-confidence SAB) (cog-confidence SNAB)))
             (ATTtv (stv ATTs ATTc)))
        (if (< 0 ATTc) (cog-merge-hi-conf-tv! ATT ATTtv)))))

(define-public (create-attr-lns TYPE)
    ;; get patient atoms and run the deduction in batch
    (cog-logger-info "Generating Attraction Links")
    ;;apply fc to get the relationship between go's and patients
    (let* ((atoms (cog-get-atoms TYPE))
            (batch-num 0)
            (batch-size (round (/ (length atoms) (current-processor-count))))
            (batch-ls (split-lst atoms batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (a)
                  (subset-attraction-rule a)) (cdr batch))) batches)
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