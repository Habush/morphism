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

(define (subset-condition-negation TYPE1 TYPE2)
  (let* ((A (Variable "$A"))
         (B (Variable "$B")))
    (Bind
      (VariableSet
        (TypedVariable A TYPE1)
        (TypedVariable B TYPE2))
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

;; Patients rule
(define (subset-negation-patients-rule)
  (subset-condition-negation (Type "PatientNode") (Type "SatisfyingSetScopeLink")))

;; Genes rule
(define (subset-condition-negation-genes-rule)
  (subset-condition-negation (Type "GeneNode") (TypeInh "ConceptNode")))