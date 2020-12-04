#!/usr/bin/racket
#lang racket
(require net/url)
(require json)
(require net/http-client)
(require racket/format)

;;; sample of how to get predicates
;;; by running /predicates you can extra a json objects with the following predicates
;;; 1.) gene_to_disease_association
;;; 2.) chemical_to_disease_or_phenotypic_feature_association
;;; 3.) diseasee_to_phenotypic_association

;;; The above pedicates math the following biolink entitites
;;; 1.) gene
;;; 2.) drug
;;; 3.) disease
;;; 4.) phenotypicfeature
(define (get-predicates url)
   (call/input-url (string->url url)
                   get-pure-port
                   (compose string->jsexpr port->string)))


;;; sample of how to build a query for chp
;;; constructs a json query object that can take in a survival time, a disease, and a single drug
;;; Input: a single drug
;;; Output: a query graph that asks this probabilistic question
;;; P(survival_time > x | drug = d1 and disease = breast cancer)
(define (build-query st disease drug)
    (define reasoner-std (make-hash))


    (define knowledge-graph (make-hash))
    (define knodes (make-hash))
    (define kedges (make-hash))
    (hash-set*! knowledge-graph (string->symbol "nodes") knodes 
                                (string->symbol "edges") kedges)

    (define node_bindings (make-hash))
    (define edge_bindings (make-hash))
    (define results (hash(string->symbol "node_bindings") node_bindings 
                         (string->symbol "edge_bindings") edge_bindings))

    (define query-graph (make-hash))
    (define qnodes (hash (string->symbol "n0") (hash
                                                (string->symbol "category")         "biolink:Gene")
                         (string->symbol "n1") (hash
                                                (string->symbol "category")         "biolink:Drug"
                                                (string->symbol "id")               drug)
                         (string->symbol "n2") (hash
                                                (string->symbol "category")         "biolink:Disease"
                                                (string->symbol "id")               disease)
                         (string->symbol "n3") (hash
                                                (string->symbol "category")         "biolink:PhenotypicFeature"
                                                (string->symbol "id")               "EFO:0000714")))

    (define qedges (hash (string->symbol "e0") (hash
                                                (string->symbol "predicate")        "biolink:GeneToDiseaseAssociation"
                                                (string->symbol "subject")          "n0"
                                                (string->symbol "object")           "n2")
                         (string->symbol "e1") (hash
                                                (string->symbol "predicate")        "biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation"
                                                (string->symbol "subject")          "n1"
                                                (string->symbol "object")           "n2")
                         (string->symbol "e2") (hash
                                                (string->symbol "predicate")        "biolink:DiseaseToPhenotypicFeatureAssociation"
                                                (string->symbol "subject")          "n2"
                                                (string->symbol "object")           "n3"
                                                (string->symbol "properties")       (hash
                                                                                        (string->symbol "qualifier")  ">="
                                                                                        (string->symbol "days")      st))))

    (hash-set! query-graph (string->symbol "nodes")      qnodes)
    (hash-set! query-graph (string->symbol "edges")      qedges)

    (hash-set*! reasoner-std (string->symbol "query_graph")     query-graph
                             (string->symbol "knowledge_graph") knowledge-graph
                             (string->symbol "results")         (list results))


    (define payload (hash))
    (hash-set payload (string->symbol "message") reasoner-std))

;;; Constructing a query and pinging CHP
;;; you can use the commneted out functionality to check which drugs are available.
;;; disease and drug are passed in as a tuple shown below. 
;;; currently only breast caner can be used as a disease
(define disease "MONDO:0007254")    ;;;Breast Cancer
(define survival-time 1000)
(define drug "CHEMBL:CHEMBL88")     ;;;CYCLOPHOSPHAMIDE

(define payload (jsexpr->bytes(build-query survival-time disease drug)))

(define conn (http-conn-open "chp.thayer.dartmouth.edu"))

(define-values (status headers in) 
    (http-conn-sendrecv! conn
               "/query/"
               #:method "POST"
               #:headers (list "Content-Type: application/json" "Accept: application/json")
               #:data payload
               #:close? #t))

(define chp-res (string->jsexpr(port->string in)))
 
;;; Extract sensitive Genes
;;; The very first result will be for the predicted survivability. 
;;; Every result thereafter contains a gene with its respective sensitivity
;;; Sensitivity values range between -1 and 1
;;; Genes closer to -1 can be thought of as having contributed more to the false assignment of P(survival_time> X | Drug^Disease)
;;; Similarly genes closer to 1 can be thought of as having contributed more to the true assignment
;;; Gene sensitivities are ordered by there absolute 
(define knowledge-graph (dict-ref (dict-ref chp-res (string->symbol "message")) (string->symbol "knowledge_graph")))
(define query-graph (dict-ref (dict-ref chp-res (string->symbol "message")) (string->symbol "query_graph")))
(define results (dict-ref (dict-ref chp-res (string->symbol "message")) (string->symbol "results")))

;;; extract probability of survival
(define survival-result (first results))


(define probability "")

(for ([qge_id (dict-keys(dict-ref survival-result (string->symbol "edge_bindings")))])
    (define kge_id "")

    (if 
        ;;; condition
        (equal? (dict-ref (dict-ref (dict-ref query-graph (string->symbol "edges")) qge_id) (string->symbol "predicate")) "biolink:DiseaseToPhenotypicFeatureAssociation")
        ;;; then
        (begin
            (set! kge_id (dict-ref(first(dict-ref (dict-ref survival-result (string->symbol "edge_bindings"))qge_id))(string->symbol "id")))
            (set! probability (dict-ref (dict-ref (dict-ref knowledge-graph (string->symbol "edges"))(string->symbol kge_id))(string->symbol "has_confidence_level")))
        )
        ;;; else
        (void)
    )
)
;;; display predicted probability     
(displayln (string-append "P(survival_time > " (~v survival-time) " | drug & disease):" (~v probability)))

;;; extract sensitive gene rankings
(define sensitivity_results (rest results))

(define genes (list))
(for ([sr sensitivity_results])
    (define gene_curie "")
    (define gene_weight "")
    (define kge_id "")
    (define sensitivity "")
    (for ([qge_id (dict-keys (dict-ref sr (string->symbol "edge_bindings")))])
        (if
            ;;; condition
            (equal? (dict-ref (dict-ref (dict-ref query-graph (string->symbol "edges")) qge_id) (string->symbol "predicate")) "biolink:GeneToDiseaseAssociation")
            ;;; then
            (begin
                (set! kge_id (dict-ref(first (dict-ref (dict-ref sr (string->symbol "edge_bindings")) qge_id)) (string->symbol "id")))
                (set! sensitivity (dict-ref (dict-ref knowledge-graph (string->symbol "edges"))(string->symbol kge_id)))
                (set! gene_curie (dict-ref sensitivity (string->symbol "subject")))
                (set! gene_weight (dict-ref sensitivity (string->symbol "value")))      
            )
            ;;; else
            (void)
        )
    )

    (define kgn_id "")
    (define gene_name "")
    (for ([qgn-id (dict-keys (dict-ref sr (string->symbol "node_bindings")))])
        (if
            ;;; condition
            (equal? (dict-ref (dict-ref (dict-ref query-graph (string->symbol "nodes")) qgn-id) (string->symbol "category")) "biolink:Gene")
            ;;; then
            (begin
                (set! kgn_id (dict-ref(first(dict-ref (dict-ref sr (string->symbol "node_bindings"))qgn-id))(string->symbol "id")))
                (set! gene_name (dict-ref(dict-ref (dict-ref knowledge-graph (string->symbol "nodes")) (string->symbol kgn_id))(string->symbol "name")))
            )
            ;;; else
            (void)
        )
    )
    (display gene_name)
    (display "\t")
    (display gene_curie)
    (display "\t")
    (display gene_weight)
    (display "\n")    
)
