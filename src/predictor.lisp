;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 20:10:56 2008 (CEST)
;; Last-Updated: 2011-08-11
;;     Update #: 13

(in-package :cltmhmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;results to .mht spec
(defun resume-zones (prediction)
  (let ((length (length prediction))
        (output (make-string 0)))
    (do ((cur-pos 0 next-pos)
         (cur-label)
         (next-pos 0))
        ((= cur-pos length) output)
      (setq cur-label (aref prediction cur-pos))
      (unless (setq next-pos (position-if-not #'(lambda (x) (char= cur-label x)) prediction :start cur-pos))
        (setq next-pos length))
      (setq output (concatenate 'string output
                                (format nil "~4T~a~18T~d~24T~d~%"
                                        (cond
                                          ((char= cur-label +membrane-labell+) "TMhelix")
                                          ((char= cur-label +cytoplasmatic-label+) "inside")
                                          ((char= cur-label +non-cytoplasmatic-label+) "output"))
                                        (1+ cur-pos) next-pos)))
      (setq cur-pos next-pos))))

;;;results to detailed spec with posterior probabilities
(defun detailed-zones (prediction thelabels post-probs)
  (let ((length (length prediction))
        (output (format nil "~37T~a~44T~a~49T~a~%"
                        (aref thelabels 0) (aref thelabels 1) (aref thelabels 2))))
    (dotimes (i length output)
      (setq output (concatenate 'string output
                                (format nil "~4T~a~36T; ~3$~45T~3$~50T~3$~%"
                                        (cond
                                          ((char= (aref prediction i) +membrane-labell+) "M")
                                          ((char= (aref prediction i) +cytoplasmatic-label+) "i")
                                          ((char= (aref prediction i) +non-cytoplasmatic-label+) "o"))
                                        (aref post-probs i 0) (aref post-probs i 1) (aref post-probs i 2)))))))

;;;predict transmembrane protein topology
(defun predict-tmpt (observation &optional identifier description)
  (let ((length (length observation))
        (obs-c (cbook *cltmhmm* observation))
        (no-helices-predicted nil))
    (multiple-value-bind (prediction)
        (viterbi-log *cltmhmm* obs-c
                     :labeled T :ending-with (list +cytoplasmatic-label+ +non-cytoplasmatic-label+))
      (setq no-helices-predicted (no-helices prediction))
      (values
       (format nil
               "# Identifier: ~a~%# Description: ~a~%# Length: ~d~%# Number of predicted TMHs: ~d~&~a~2%"
               (if identifier identifier "") (if description description "") length no-helices-predicted
               (resume-zones prediction))
       prediction
       (multiple-value-bind (thelabels post-probs)
           (posterior-probs-labels *cltmhmm* obs-c '(+cytoplasmatic-label+ +non-cytoplasmatic-label+))
           (detailed-zones prediction thelabels post-probs))))))

