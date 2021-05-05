(use-modules (ice-9 threads))
(use-modules (srfi srfi-1))

;; Subset rule for Genes
;;
;; Member <STV>
;;   A 
;;   B 
;; |-
;; Subset <TV>
;;   (Set A)
;;   B
;;

;; Rule for Subset
;;
(define (subset-genes-rule A)
  (let* ((B (Variable "$B")))
    (Bind
      (VariableSet
        (TypedVariable B (TypeInh "ConceptNode")))
      (Present
        (Member A B))
      (ExecutionOutput
        (GroundedSchema "scm: generate-subset-from-member")
        (List
          ;; conclusion
          (Subset (Set A) B)
          (Member A B))))))

(define (generate-subset-from-member conclusion . premises)
  (if (= (length premises) 1)
    (let* ((Ss conclusion)
          (memb (car premises))
          (st (cog-mean memb))
          (conf (cog-confidence memb)))
    (if (> conf 0) (cog-merge-hi-conf-tv! Ss (stv st conf))))))

;; Name
(define subset-genes-rule-name
  (DefinedSchemaNode "subset-genes-rule"))
(DefineLink subset-genes-rule-name subset-genes-rule)


(define-public (create-subset-lns TYPE)
    (cog-logger-info "Generating Subset Links")
    (let* ((atoms (cog-get-atoms TYPE))
            (batch-num 0)
            (batch-size (/ (len genes) (current-processor-count)))
            (batch-ls (split-lst genes batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls)))
        
        (n-par-for-each (current-processor-count)  (lambda (batch)
              (for-each (lambda (a)
                  (subset-genes-rule a)) batch)) batches)
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