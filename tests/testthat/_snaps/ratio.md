# ratio.estimates

    Code
      ratio(y_num = apisrs$api.stu, y_den = apisrs$enroll, xsample = apisrs$stype,
      xpop = apipop$stype, pi = apisrs$pw^(-1), estimator = "postStrat", var_est = TRUE,
      var_method = "LinHB", datatype = "raw")
    Output
      $ratio_est
      [1] 0.8258652
      
      $ratio_var_est
      [1] 8.851525e-05
      

