# XCut

This is XCut, an expander decomposition based solver for normalized cut.

## Building

The project is built using [Bazel](https://bazel.build). build system.

``` shell
# Build optimized configuration
bazel build -c opt //main:rw

# Build debug configuration
bazel build -c dbg //main:rw

# Build with better debug information and sanitizers
bazel build --config debug --compilation_mode dbg --sandbox_debug //main:rw
```

## Running

Invoke `rw`, the options are documented in the file `main/rw.cpp`. For computing a
normalized cut on a graph, consider the following example invocation:

``` shell
./bazel-bin/main/rw -logtostderr -hierarchy -experimentType 3 -numParts 30 graph-file.mtx
```

For `chaco` files, use the `-chaco` flag.