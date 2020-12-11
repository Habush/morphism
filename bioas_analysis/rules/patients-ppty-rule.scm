;; Patients Infered ppty Subset introduction rule 
;;
;; SubsetLink
;;   And
;;     Go term or Pathway
;;     profiled-genes
;;   SatisfyingSetScopeLink
;;     VariableNode "$G"
;;     EvaluationLink
;;       LazyExecutionOutputLink
;;         SchemaNode "make-overexpression-predicate-for-gene"
;;         VariableNode "$G"
;;       Patient
;; |-
;; Subset 
;;   Set Patient
;;   SatisfyingSetScopeLink
;;     Variable "$P"
;;     SubsetLink
;;       And
;;         Go term or Pathway
;;         profiled-genes
;;       SatisfyingSetScopeLink
;;         VariableNode "$G"
;;         EvaluationLink
;;           LazyExecutionOutputLink
;;             SchemaNode "make-overexpression-predicate-for-gene"
;;             VariableNode "$G"
;;           Variable "$P"

(define patient-ppty-subset-rule
  (let* ((B (Variable "$B"))
         (P (Variable "$P"))
         (S (Variable "$S")))
    (Bind
      (VariableSet
        (TypedVariable B (TypeInh "ConceptNode"))
        (TypedVariable P (Type "PatientNode"))
        (TypedVariable S (Type "SchemaNode")))
      (Present
        (SubsetLink
            (And
                B
                (ConceptNode "profiled-genes"))
            (SatisfyingSetScopeLink
                (VariableNode "$G")
                (EvaluationLink
                (LazyExecutionOutputLink
                    S
                    (VariableNode "$G"))
                P))))
      (ExecutionOutput
        (GroundedSchema "scm: generate-subset")
        (List
          ;; conclusion
          (Subset
            (Set P)
            (SatisfyingSetScopeLink
              (Variable "$Pt")
              (SubsetLink
                (And
                  B
                  (ConceptNode "profiled-genes"))
                (SatisfyingSetScopeLink
                  (VariableNode "$G")
                  (EvaluationLink
                    (LazyExecutionOutputLink
                      S
                      (VariableNode "$G"))
                    (Variable "$Pt"))))))
          ;; premises
          (SubsetLink
            (And
                B
                (ConceptNode "profiled-genes"))
            (SatisfyingSetScopeLink
                (VariableNode "$G")
                (EvaluationLink
                (LazyExecutionOutputLink
                    S
                    (VariableNode "$G"))
                P))))))))

(define (generate-subset conclusion . premises)
;(ure-logger-debug "(conclusion=~a . premises=~a)" conclusion premises)
  (if (= (length premises) 1)
    (let* ((Ss conclusion)
          (pr (car premises))
          (st (cog-mean pr))
          (conf (cog-confidence pr)))
      (if (> conf 0) (cog-merge-hi-conf-tv! Ss (stv st conf))))))

; Name the rule
(define patient-ppty-subset-rule-name
  (DefinedSchemaNode "patient-ppty-subset-rule"))
(DefineLink patient-ppty-subset-rule-name
  patient-ppty-subset-rule)
