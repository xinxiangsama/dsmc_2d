#include "Output.h"

Output::Output(Run *run)
{
    m_run = run;
}

void Output::Write2HDF5(const std::string &filename)
{
    auto N1local = m_run->m_mesh->getnumberCellsX();
    auto N2local = m_run->m_mesh->getnumberCellsY();
    auto N1global = m_run->m_mesh->getnumberCellsXGlobal();
    auto N2global = m_run->m_mesh->getnumberCellsYGlobal();
    auto offsetX = m_run->m_mesh->getoffsetX();
    auto offsetY = m_run->m_mesh->getoffsetY();

    H5::FileAccPropList fapl;
    H5Pset_fapl_mpio(fapl.getId(), MPI_COMM_WORLD, MPI_INFO_NULL);

    H5::H5File file(filename, H5F_ACC_TRUNC, H5::FileCreatPropList::DEFAULT, fapl);
    hsize_t global_dims[2] = {N1global, N2global};
    hsize_t local_dims[2] = {N1local, N2local};
    hsize_t offset[2] = {offsetX, offsetY};

    H5::DataSpace global_dataspace(2, global_dims);
    H5::DataSpace local_dataspace(2, local_dims);
    global_dataspace.selectHyperslab(H5S_SELECT_SET, local_dims, offset);

    H5::DSetCreatPropList plist;
    plist.setChunk(2, local_dims);

    H5::DataSet dataset_U = file.createDataSet("U", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_V = file.createDataSet("V", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_P = file.createDataSet("P", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_T = file.createDataSet("T", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_Rho = file.createDataSet("Rho", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);

    std::vector<double> U(N1local * N2local, 0.0);
    std::vector<double> V(N1local * N2local, 0.0);
    std::vector<double> P(N1local * N2local, 0.0);
    std::vector<double> T(N1local * N2local, 0.0);
    std::vector<double> Rho(N1local * N2local, 0.0);
    for (size_t i = 0; i < N1local; ++i){
        for (size_t j = 0; j < N2local; ++j){
            auto cell = m_run->m_cells[i * N2local + j].get();
            auto phase = cell->getphase();
            auto element = cell->getelement();

            U[i * N2local + j] = phase->getvelocity().x();
            V[i * N2local + j] = phase->getvelocity().y();
            P[i * N2local + j] = phase->getpressure();
            T[i * N2local + j] = phase->gettemperature();
            Rho[i * N2local + j] = phase->getdensity();
        }
    }

    dataset_U.write(U.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_V.write(V.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_P.write(P.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_T.write(T.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_Rho.write(Rho.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
}
