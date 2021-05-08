(use-modules (ice-9 threads))
(use-modules (srfi srfi-1))

;; Subset Negation rule
;;
;; Subset <STV>
;;   Set A <ATV> 
;;   B <BTV>
;; |-
;; Subset <TV>
;;   Not (Set A)
;;   B
;;

(define (subset-condition-negation A)
  (let* ((B (Variable "$B")))
    (Bind
      (VariableSet
        (TypedVariable B (TypeInh "ConceptNode")))
      (Present
        (Subset (Set A) B))
      (ExecutionOutput
        (GroundedSchema "scm: get-subset-condition-negation")
        (List
          ;; Conclusion
          (Subset (Not (Set A)) B)
          ;; Premises
          (Subset (Set A) B)
          (Set A)
          B)))))

;; Formula
(define (get-subset-condition-negation conclusion . premises)
 ;(ure-logger-debug "(conclusion=~a . premises=~a)" conclusion premises)
  (if (= (length premises) 3)
      (let* ((NS conclusion)
             (S (car premises))
             (A (cadr premises))
             (B (caddr premises))
             (Ss (cog-mean S))
	           (Sc (cog-confidence S))
             (As (cog-mean A))
             (Ac (cog-confidence A))
             (Bs (cog-mean B))
	           (NAs (- 1 As))
             (NSs (if (< As 1)
                      (/ (- Bs (* Ss As)) NAs)
                      1))
             (NSc (if (< As 1)
                      (min (count->confidence (* (confidence->count Ac) NAs)) Sc)
                      0))
             (NStv (stv NSs NSc)))
        (cog-merge-hi-conf-tv! NS NStv))))

(define-public (create-subset-neg-lns TYPE)
    (cog-logger-info "Generating Subset Negation Links")
    (let* ((atoms (cog-get-atoms TYPE))
            (batch-num 0)
            (batch-size (round (/ (length atoms) (current-processor-count))))
            (batch-ls (split-lst atoms batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (a)
                  (subset-condition-negation a)) (cdr batch))) batches)
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