;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 20:18:05 2008 (CEST)
;; Last-Updated: 2011-08-11
;;     Update #: 14

(in-package :cltmhmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-next-helix (sequence &optional (start 0) (tm-label #\M))
  (let* ((first-tm (position tm-label sequence :start start))
         (last-tm (when first-tm (position-not tm-label sequence :start first-tm))))
    (if (not first-tm)
        (return-from find-next-helix (values 0))
        (progn
          (unless last-tm
            (setq last-tm (length sequence)))
          (return-from find-next-helix
            (values (- last-tm first-tm) first-tm last-tm))))))

(defun no-helices (sequence &optional (start 0) (tm-label #\M))
  (let ((leng 0) (end start) (no 0))
    (while t
      (multiple-value-setq (leng start end) (find-next-helix sequence end tm-label))
      (cond
        ((= leng 0) (return no))
        (t (incf no))))))

(defun tm-prediction-accuracy (original prediction &key (min-coincidence 5) (tm-label #\M) (verbose nil))
  (let* ((size (length original)) ;of course, .prediction. is supposed to have same length
         (orig-ort (aref original 0))
         (pred-ort (aref prediction 0))
         (correct-ort (equal orig-ort pred-ort))
         (no-helices 0)
         (no-helices-pred 0)
         (no-helices-shifted 0)
         (positives 0)
         (false-positives 0)
         (false-negatives 0)
         (false-splits 0)
         (false-merges 0))
     ;;transverse the original
    (let ((leng 0) (start 0) (end 0))
      (while t
        (multiple-value-setq (leng start end) (find-next-helix original end tm-label))
        (cond
          ((= 0 leng) (return))
          (t (incf no-helices)
             (setq leng (count tm-label prediction :start start :end end))
             (cond
               ((= leng 0) (incf false-negatives))
               (t nil))
             (setq leng (no-helices (subseq prediction start end) 0 tm-label)) ;let me use this var
             (incf false-splits (- (if (zerop leng) 1 leng)  1))))))
    ;;transverse the prediction
    (let ((leng 0) (start 0) (end 0))
      (while t
        (multiple-value-setq (leng start end) (find-next-helix prediction end tm-label))
        (cond
          ((= 0 leng) (return))
          (t (incf no-helices-pred)
             (setq leng (count tm-label original :start start :end end))
             (cond
               ((>= leng min-coincidence) (incf positives))
               ((and (>= leng 1) (= 1 (setq leng (no-helices (subseq original start end) 0 tm-label))))
                (incf no-helices-shifted))
               ((> leng 1) (incf false-merges (- leng 1)))
               (t (incf false-positives)))))))
    (when verbose
      (format t "length: ~a~%correct orientation: ~a (orig: ~a, pred: ~a)~%no helices: ~a~%no helices predicted: ~a~%   correct: ~a~%   false positives: ~a~%   false negatives: ~a~%   false splits: ~a~%   false merges: ~a~%   helices shifted: ~a~%"
              size correct-ort orig-ort pred-ort no-helices no-helices-pred positives false-positives false-negatives false-splits false-merges no-helices-shifted))
    (list size orig-ort pred-ort correct-ort no-helices no-helices-pred
          positives false-positives false-negatives false-splits false-merges no-helices-shifted)))

;;;test the accuracy for multiple predictions
(defun multiple-tm-prediction-accuracy (originals predictions
                                        &key (min-coincidence 5) (tm-label #\M) (verbose nil))
  (macrolet ((an (e &optional (list 'output))
               `(car (push ,e ,list))))
    (let ((proteins 0)
          (sp-proteins 0)
          (correct-top 0)
          (inverted-top 0)
          (correct-n-terminal 0)
          (over-predictions 0)
          (sp-over-predictions 0) ;single spanning
          (under-predictions 0)
          (sp-under-predictions 0)
          (over&under-predictions 0)
          (sp-over&under-predictions 0)
          (no-real-helices 0)
          (false-positives 0)
          (false-negatives 0)
          (shifted-helices 0)
          (merged-helices 0)
          (splited-helices 0)
          (output 0))
      (loop
         for o in originals
         for p in predictions
         for i = 0 then (1+ i) do
         ;;length, original orientation, prediction orientation, correct orientation, real helices, predicted helices
         ;;correct helices, false positives, false negatives, false splits, false merges, helices shited
           (destructuring-bind (l oo po c rh ph ch fp fn fs fm sh) (tm-prediction-accuracy o p
                                                                                         :min-coincidence min-coincidence
                                                                                         :tm-label tm-label
                                                                                         :verbose verbose)
             (declare (ignore l oo po))
             (when verbose
               (format t "~a~%~a~3%" o p))
             (cond
               ((and c (= rh ph ch)) (incf correct-top))
               ((= rh ph ch) (incf inverted-top))
               ((> fp 0) (cond
                           ((> fn 0) (incf over&under-predictions) (when (= 1 rh) (incf sp-over&under-predictions)))
                           (t (incf over-predictions) (when (= 1 rh) (incf sp-over-predictions)))))
               ((> fn 0) (incf under-predictions) (when (= 1 rh) (incf sp-under-predictions)))
               (t))
             (when c (incf correct-n-terminal))
             (when (= rh 1) (incf sp-proteins))
             (incf no-real-helices rh)
             (incf false-positives fp)
             (incf false-negatives fn)
             (incf shifted-helices sh)
             (incf merged-helices fm)
             (incf splited-helices fs))
         finally
           (setf proteins (1+ i))
           (print over-predictions)
           (format t "
~{Number of proteins: ~50T~d
~5Tof which single-spanning: ~50T~d ~56T~d%
Correctly predicted topology: ~50T~d ~56T~d%
~5TInvertedly predicted topology: ~50T~d ~56T~d%
Correctly predicted N-terminal: ~50T~d ~56T~d%
Under-predictions: ~50T~d ~56T~d%
~5Tof which single-spanning: ~50T~d ~56T~d%
Over-predictions: ~50T~d ~56T~d%
~5Tof which single-spanning: ~50T~d ~56T~d%
Both over- and under-predictions: ~50T~d ~56T~d%
~5Tof which single-spanning: ~50T~d ~56T~d%

Total number of real helices: ~50T~d
Number of over-predicted helices: ~50T~d ~56T~d%
Number of under-predicted helices: ~50T~d ~56T~d%
Number of shifted helix predictions: ~50T~d ~56T~d%
Number of falsely merged helices: ~50T~d ~56T~d%
Number of falsely split helices: ~50T~d ~56T~d% ~} ~%"
               (list (an proteins)
                     (an sp-proteins) (%100 sp-proteins proteins)
                     (an correct-top) (%100 correct-top proteins)
                     (an inverted-top) (%100 inverted-top proteins)
                     (an correct-n-terminal) (%100 correct-n-terminal proteins)
                     (an under-predictions) (%100 under-predictions proteins)
                     (an sp-under-predictions) (%100 sp-under-predictions proteins)
                     (an over-predictions) (%100 over-predictions proteins)
                     (an sp-over-predictions) (%100 sp-over-predictions proteins)
                     (an over&under-predictions) (%100 over&under-predictions proteins)
                     (an sp-over&under-predictions) (%100 sp-over&under-predictions proteins)

                     (an no-real-helices)
                     (an false-positives) (an (%100 false-positives no-real-helices))
                     (an false-negatives) (an (%100 false-negatives no-real-helices))
                     (an shifted-helices) (an (%100 shifted-helices no-real-helices))
                     (an merged-helices) (an (%100 merged-helices no-real-helices))
                     (an splited-helices) (an (%100 splited-helices no-real-helices))))
       (return output)))))


