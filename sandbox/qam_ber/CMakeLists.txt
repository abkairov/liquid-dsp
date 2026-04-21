add_executable(qam_ber EXCLUDE_FROM_ALL
    qam_ber.c
    modem_pipeline.c
    channel_sim.c)
target_include_directories(qam_ber PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(qam_ber liquid)
