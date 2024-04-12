load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")


cc_library(
  name = "input_util",
  hdrs = ["util.hpp"],
  deps = [
    "//lib:cluster_util",
  ],
)

cc_binary(
  name = "rw",
  srcs = ["rw.cpp"],
  deps = [
    "input_util",
    "//lib:cluster_util",
    "@com_google_glog//:glog",
  ]
)

cc_binary(
  name = "flow",
  srcs = ["flow.cpp"],
  deps = [
    "input_util",
    "//lib:cluster_util",
    "@com_google_glog//:glog",
  ]
)

cc_binary(
  name = "quality",
  srcs = ["quality.cpp"],
  deps = [
    "input_util",
    "//lib:cluster_util",
    "@com_google_glog//:glog",
  ]
)

cc_binary(
  name = "eval_ncut",
  srcs = ["eval_ncut.cpp"],
  deps = [
    "input_util",
    "//lib:cluster_util",
    "@com_google_glog//:glog",
  ]
)