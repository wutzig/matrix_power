from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps

class snake_tracer(ConanFile):
    name = "matrix_power"
    version = 0.1
    description = "profiliing BLAS algorithms"
    author = "Michael"
    settings = "os", "arch", "compiler", "build_type"

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.21]")
    
    def requirements(self):
        self.requires("openblas/[>=0.3]")

    def layout(self):
        self.folders.source = "src"
        self.folders.build = "build"
        self.folders.generators = "cmake"

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

        # creating .cmake files for find_package dependencies
        deps = CMakeDeps(self)
        deps.generate()
    
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
            