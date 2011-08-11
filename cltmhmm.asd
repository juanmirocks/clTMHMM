;; Author: Juan Miguel Cejuela
;; Created: Mon Sep 22 11:02:11 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 8

(defpackage :net.ashrentum.cltmhmm-system (:use :asdf :cl))
(in-package :net.ashrentum.cltmhmm-system)

(defsystem cltmhmm
    :name "cltmhmm"
    :author "Juan Miguel Cejuela"
    :version "0.2"
    :maintainer "Juan Miguel Cejuela"
    :licence "GPL 3"
    :description "TMHMM-based TM Protein Topology Predictor in CL"
    :components
    ((:module
      "src"
      :components
      ((:file "packages")

       (:file "model"
        :depends-on ("packages"))

       (:file "accuracy-validation"
        :depends-on ("packages"))

       (:file "predictor"
        :depends-on ("packages" "accuracy-validation" "model"))

       (:file "other-results-interpreter"
        :depends-on ("packages")))))

    :depends-on
    (:jmc.cl.utils
     :cl-hmm))

