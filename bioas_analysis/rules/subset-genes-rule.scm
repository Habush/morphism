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
(define subset-genes-rule
  (let* ((A (Variable "$A"))
         (B (Variable "$B")))
    (Bind
      (VariableSet
        (TypedVariable A (Type "GeneNode"))
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
