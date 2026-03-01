cmake_minimum_required(VERSION 3.10)
project(CurveIntersection)

# 1. 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------
# 2. 手动定义库路径 (请根据你电脑上的实际路径修改)
# ---------------------------------------------------------

# Eigen 路径 (只需要头文件)
set(EIGEN_INCLUDE_DIR "D:/lib/eigen")

# OpenGJK 路径 (需要头文件 + 库文件)
set(OPENGJK_INCLUDE_DIR "D:/lib/openGJK/include")
set(OPENGJK_LIB_DIR     "D:/lib/openGJK/build/Debug")
set(OPENGJK_LIBRARY     "opengjk_ce.lib") 
set(OPENGJK_DLL_PATH    "${OPENGJK_LIB_DIR}/opengjk_ce.dll")

# ---------------------------------------------------------
# 3. 查找源文件并生成目标
# ---------------------------------------------------------
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# ---------------------------------------------------------
# 4. 包含头文件目录 (Include)
# ---------------------------------------------------------
target_include_directories(${PROJECT_NAME} PRIVATE 
    ${EIGEN_INCLUDE_DIR}
    ${OPENGJK_INCLUDE_DIR}
)

# ---------------------------------------------------------
# 5. 链接库文件 (Link)
# ---------------------------------------------------------
target_link_directories(${PROJECT_NAME} PRIVATE ${OPENGJK_LIB_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGJK_LIBRARY})

# ---------------------------------------------------------
# 6. 核心：将 Example 文件夹转移到生成的 Release/Debug 目录
# ---------------------------------------------------------
# 使用 POST_BUILD 确保每次编译后都会同步最新的 Example 文件
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/Example"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${OPENGJK_DLL_PATH}"
    "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
    COMMENT "Moving Example files to: $<TARGET_FILE_DIR:${PROJECT_NAME}>/Example"
)