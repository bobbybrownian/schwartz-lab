class NewstuffsController < ApplicationController

  def new
    @newstuff = Newstuff.new
    @newstuffs = Newstuff.all
  end

  def create

    # grab params for newstuff
    newstuff_params = params["newstuff"]

    description = newstuff_params["description"]
    @newstuff = Newstuff.new(
                          description: description)
    @newstuff.save
    # binding.pry
    redirect_to maintenance_index_path
  end

  def edit

  end

  def updated

  end

  def destroy
    newstuff = newstuff.find(params[:id])
    newstuff.destroy
    redirect_to new_newstuff_path
  end

  end
