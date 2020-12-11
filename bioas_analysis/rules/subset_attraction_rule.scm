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

(define (subset-attraction-rule TYPE1 TYPE2)
  (let* ((A (Variable "$A"))
         (B (Variable "$B")))
    (BindLink
      (VariableSet
        (TypedVariable A TYPE1)
        (TypedVariable B TYPE2))
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

(define subset-attraction-patients-rule
  (subset-attraction-rule (Type "PatientNode") (Type "SatisfyingSetScopeLink")))
(define subset-attraction-patients-rule-name
  (DefinedSchemaNode "subset-attraction-patients-rule"))
(DefineLink subset-attraction-patients-rule-name
  subset-attraction-patients-rule)

;; Genes
(define subset-attraction-genes-rule
  (subset-attraction-rule (Type "GeneNode") (TypeInh "ConceptNode")))
; Name the rule
(define subset-attraction-genes-rule-name
  (DefinedSchemaNode "subset-attraction-genes-rule"))
(DefineLink subset-attraction-genes-rule-name
  subset-attraction-genes-rule)