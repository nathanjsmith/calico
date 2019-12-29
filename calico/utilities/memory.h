// Copyright (c) 2019, Nathan Smith <nathanjsmith@gmail.com>
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#ifndef __CALICO__UTILITIES__MEMORY_H__
#define __CALICO__UTILITIES__MEMORY_H__

#include <cstddef>

// Windows
#ifdef _WIN32
#include <malloc.h>

namespace calico {
namespace utilities {
namespace memory {

void* aligned_alloc(std::size_t alignment, std::size_t size) {
    return _aligned_malloc(size, alignment);
}

void aligned_free(void *memory) {
    _aligned_free(memory);
}

}// end namespace calico::utilities::memory
}// end namespace calico::utilities
}// end namespace calico

#else
// Linux
#include <stdlib.h>

namespace calico {
namespace utilities {
namespace memory {

void* aligned_alloc(std::size_t alignment, std::size_t size) {
    return ::aligned_alloc(size, alignment);
}

void aligned_free(void *memory) {
    ::aligned_free(memory);
}

}// end namespace calico::utilities::memory
}// end namespace calico::utilities
}// end namespace calico
#endif

#endif // end include guard