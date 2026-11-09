Allocation failure regression tests can be enabled with:

```sh
cmake -S /path/to/RingDecomposerLib -B /tmp/rdl-tests -DRDL_BUILD_ALLOCATION_TESTS=ON
cmake --build /tmp/rdl-tests
ctest --test-dir /tmp/rdl-tests --output-on-failure
```

The test library substitutes allocation functions without changing the production
library. Each operation is repeated with successive allocations forced to fail.
Tests check error returns, diagnostic callbacks, output pointers, outstanding
allocations, and retrying edge insertion after failure. Fixtures cover dense
graphs, disconnected components, even cycles, long paths, and empty graphs.

For AddressSanitizer and UndefinedBehaviorSanitizer, configure a separate debug
build with `-DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"`.
The test allocator also detects leaks on platforms without LeakSanitizer.
