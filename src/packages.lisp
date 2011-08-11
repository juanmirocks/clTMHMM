;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 20:22:28 2008 (CEST)
;; Last-Updated: 2011-08-11
;;     Update #: 11

(defpackage :net.ashrentum.cltmhmm
  (:use :cl :jmcejuela :cl-hmm)
  (:nicknames :cltmhmm)
  (:export

   :tm-prediction-accuracy
   :multiple-tm-prediction-accuracy

   :predict-tmpt))

