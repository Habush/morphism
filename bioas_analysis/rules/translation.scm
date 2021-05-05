(define-module (pln-bio rules translation)
   #:use-module (opencog) 
   #:use-module (opencog exec) 
   #:use-module (opencog ure) 
   #:use-module (opencog logger)
   #:use-module (srfi srfi-1)
   #:use-module (pln-bio rules rule-utils)
)
;; Crisp rules about translating a link into another link

;; Helpers
(define ConceptT (TypeInh "ConceptNode"))
(define GeneT (Type "GeneNode"))

;; Inheritance to Subset
(define present-inheritance-to-subset-translation-rule
  (gen-present-link-translation-rule InheritanceLink SubsetLink ConceptT))
(define present-inheritance-to-subset-translation-rule-name
  (DefinedSchemaNode "present-inheritance-to-subset-translation-rule"))
(DefineLink present-inheritance-to-subset-translation-rule-name
  present-inheritance-to-subset-translation-rule)
  
